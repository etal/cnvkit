"""Segmentation by Hidden Markov Model.

Pure numpy/scipy implementation with distance-dependent transitions and
joint (log2, BAF) emissions parameterized by (purity, ploidy).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.special import betaln, gammaln, logsumexp

from ..cnary import CopyNumArray as CNA
from ..descriptives import biweight_midvariance
from ..segfilters import squash_by_groups

if TYPE_CHECKING:
    from cnvlib.vary import VariantArray
    from numpy.typing import NDArray

# Configuration constants
DEFAULT_LENGTH_SCALE = 5e6  # bp; p_stay(d) = exp(-d/L)
TRANSITION_PENALTY = 8.0  # penalty for CN-distant transitions
LOG2_STDEV_FLOOR = 0.05  # minimum stdev for log2 Gaussian
BETABINOM_RHO = 200.0  # beta-binomial concentration; higher → less overdispersion
MIN_VARIANTS_THRESHOLD = 50  # for variants_in_segment


# --- State space ---


def enumerate_states(
    max_copies: int, include_baf: bool
) -> list[tuple[int, int | None]]:
    """Enumerate HMM states as (copy_number, minor_allele) tuples.

    Parameters
    ----------
    max_copies : int
        Maximum total copy number.
    include_baf : bool
        If True, enumerate minor allele counts (0..cn//2) for each CN,
        giving a richer state space for joint log2+BAF modeling.

    Returns
    -------
    list of (int, int | None)
        Each entry is (cn, minor) where minor is None when BAF is not used.
    """
    states: list[tuple[int, int | None]] = []
    if include_baf:
        for cn in range(max_copies + 1):
            for minor in range(cn // 2 + 1):
                states.append((cn, minor))
    else:
        for cn in range(max_copies + 1):
            states.append((cn, None))
    return states


# --- Emission model ---


def _expected_log2(
    cn: NDArray[np.float64],
    purity: float,
    ploidy: float,
) -> NDArray[np.float64]:
    """Expected log2 ratio for given CN, purity, ploidy.

    Formula: log2((p*cn + (1-p)*2) / (p*ploidy + (1-p)*2))
    """
    p = purity
    numerator = p * cn + (1 - p) * 2.0
    denominator = p * ploidy + (1 - p) * 2.0
    # Floor numerator to avoid log2(0) for cn=0
    ratio = np.maximum(numerator, 1e-4) / denominator
    return np.log2(ratio)  # type: ignore[no-any-return]


def _expected_minor_freq(
    minor: NDArray[np.float64],
    cn: NDArray[np.float64],
    purity: float,
) -> NDArray[np.float64]:
    """Expected minor allele frequency (<= 0.5) for given minor allele count, CN, purity.

    Formula: (p*minor + (1-p)*1) / (p*cn + (1-p)*2)
    For cn=0, return 0.5 (no allelic information).
    """
    p = purity
    numerator = p * minor + (1 - p) * 1.0
    denominator = p * cn + (1 - p) * 2.0
    # Avoid division by zero for cn=0 at high purity
    with np.errstate(divide="ignore", invalid="ignore"):
        baf = np.where(denominator > 0, numerator / denominator, 0.5)
    return baf


def _betabinom_logpmf(
    k: NDArray[np.float64],
    n: NDArray[np.float64],
    a: NDArray[np.float64],
    b: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Log PMF of the beta-binomial distribution, vectorized.

    Uses scipy.special gammaln/betaln for numerical stability.
    When n=0, returns 0.0 (no observation contributes nothing).
    """
    return (  # type: ignore[no-any-return]
        gammaln(n + 1)
        - gammaln(k + 1)
        - gammaln(n - k + 1)
        + betaln(k + a, n - k + b)
        - betaln(a, b)
    )


def log_emission_probs(
    log2_obs: NDArray[np.float64],
    minor_counts: NDArray[np.float64] | None,
    depths: NDArray[np.float64] | None,
    states: list[tuple[int, int | None]],
    purity: float,
    ploidy: float,
    log2_stdev: float,
) -> NDArray[np.float64]:
    """Compute log emission probabilities for all observations and states.

    Parameters
    ----------
    log2_obs : ndarray (T,)
        Observed log2 ratios.
    minor_counts : ndarray (T,) or None
        Aggregated minor allele counts per bin.
    depths : ndarray (T,) or None
        Aggregated read depths per bin.
    states : list of (cn, minor|None)
        HMM state definitions.
    purity, ploidy : float
        Tumor purity and ploidy.
    log2_stdev : float
        Standard deviation for log2 Gaussian emission.

    Returns
    -------
    ndarray (T, N)
        Log probability of each observation under each state.
    """
    N = len(states)

    cn_arr = np.array([s[0] for s in states], dtype=np.float64)
    expected_l2 = _expected_log2(cn_arr, purity, ploidy)  # (N,)

    # Log2 component: N(obs | expected, stdev)
    # log N(x|mu,s) = -0.5*((x-mu)/s)^2 - log(s) - 0.5*log(2*pi)
    diff_l2 = log2_obs[:, np.newaxis] - expected_l2[np.newaxis, :]  # (T, N)
    log_p = -0.5 * (diff_l2 / log2_stdev) ** 2 - np.log(
        log2_stdev * np.sqrt(2.0 * np.pi)
    )

    # BAF component: beta-binomial on aggregated allele counts
    if minor_counts is not None and depths is not None:
        minor_arr = np.array(
            [s[1] if s[1] is not None else 0 for s in states], dtype=np.float64
        )
        expected_freq = _expected_minor_freq(minor_arr, cn_arr, purity)  # (N,)

        # Beta-binomial parameters with floor to keep betaln defined
        alpha = np.maximum(expected_freq * BETABINOM_RHO, 0.5)  # (N,)
        beta = np.maximum((1.0 - expected_freq) * BETABINOM_RHO, 0.5)  # (N,)

        k = minor_counts[:, np.newaxis]  # (T, 1)
        n = depths[:, np.newaxis]  # (T, 1)
        a = alpha[np.newaxis, :]  # (1, N)
        b = beta[np.newaxis, :]  # (1, N)

        log_baf = _betabinom_logpmf(k, n, a, b)  # (T, N)
        # Zero-depth bins contribute nothing
        zero_depth = depths == 0
        log_baf[zero_depth, :] = 0.0
        log_p += log_baf

    return log_p  # type: ignore[no-any-return]


# --- Transition model ---


def _switch_weights(
    n_states: int,
    cn_values: NDArray[np.int_],
    penalty: float,
) -> NDArray[np.float64]:
    """Precompute normalized off-diagonal switch weights.

    Weight for switching from state i to state j is proportional to
    exp(-penalty * |cn_i - cn_j|), normalized per row (excluding diagonal).

    Returns
    -------
    ndarray (N, N)
        Normalized switch weight matrix (rows sum to 1, diagonal is 0).
    """
    cn = cn_values.astype(np.float64)
    diff = np.abs(cn[:, np.newaxis] - cn[np.newaxis, :])
    weights = np.exp(-penalty * diff)
    # Zero diagonal
    np.fill_diagonal(weights, 0.0)
    # Normalize rows
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-300)
    weights /= row_sums
    return weights  # type: ignore[no-any-return]


def log_transition_matrices(
    distances: NDArray[np.float64],
    n_states: int,
    cn_values: NDArray[np.int_],
    length_scale: float = DEFAULT_LENGTH_SCALE,
    penalty: float = TRANSITION_PENALTY,
) -> NDArray[np.float64]:
    """Compute per-step log transition matrices from inter-bin distances.

    Parameters
    ----------
    distances : ndarray (T-1,)
        Inter-bin distances in bp.
    n_states : int
        Number of HMM states.
    cn_values : ndarray (N,)
        Copy number for each state.
    length_scale : float
        Distance scale for stay probability: p_stay = exp(-d/L).
    penalty : float
        Penalty for CN-distant transitions.

    Returns
    -------
    ndarray (T-1, N, N)
        Log transition probability matrices.
    """
    sw = _switch_weights(n_states, cn_values, penalty)  # (N, N)
    p_stay = np.exp(-distances / length_scale)  # (T-1,)
    # Clip to avoid log(0)
    p_stay = np.clip(p_stay, 1e-300, 1.0 - 1e-10)

    # T[t, i, j] = p_stay[t]*I(i==j) + (1-p_stay[t])*sw[i,j]
    eye = np.eye(n_states)
    # Shape: (T-1, N, N)
    trans = (
        p_stay[:, np.newaxis, np.newaxis] * eye[np.newaxis, :, :]
        + (1.0 - p_stay[:, np.newaxis, np.newaxis]) * sw[np.newaxis, :, :]
    )

    # Ensure valid probabilities
    trans = np.maximum(trans, 1e-300)
    return np.log(trans)  # type: ignore[no-any-return]


# --- Forward algorithm ---


def forward_log(
    log_emit: NDArray[np.float64],
    log_trans: NDArray[np.float64],
    log_start: NDArray[np.float64],
) -> tuple[NDArray[np.float64], float]:
    """Forward algorithm in log space.

    Parameters
    ----------
    log_emit : ndarray (T, N)
        Log emission probabilities.
    log_trans : ndarray (T-1, N, N)
        Log transition matrices.
    log_start : ndarray (N,)
        Log start probabilities.

    Returns
    -------
    alpha : ndarray (T, N)
        Forward log probabilities.
    log_likelihood : float
        Total log-likelihood of the observation sequence.
    """
    T, N = log_emit.shape
    alpha = np.empty((T, N), dtype=np.float64)
    alpha[0] = log_start + log_emit[0]

    for t in range(1, T):
        # alpha[t, j] = log_emit[t, j] + logsumexp(alpha[t-1] + log_trans[t-1, :, j])
        # Vectorize over j: alpha[t-1, :] broadcast with log_trans[t-1, :, :]
        # log_trans[t-1] is (N, N) where [i, j] = log P(i->j)
        # We want for each j: logsumexp_i(alpha[t-1, i] + log_trans[t-1, i, j])
        alpha[t] = log_emit[t] + logsumexp(
            alpha[t - 1, :, np.newaxis] + log_trans[t - 1], axis=0
        )

    log_likelihood = float(logsumexp(alpha[-1]))
    return alpha, log_likelihood


# --- Viterbi algorithm ---


def viterbi_log(
    log_emit: NDArray[np.float64],
    log_trans: NDArray[np.float64],
    log_start: NDArray[np.float64],
) -> NDArray[np.int_]:
    """Viterbi algorithm with backtracking in log space.

    Parameters
    ----------
    log_emit : ndarray (T, N)
        Log emission probabilities.
    log_trans : ndarray (T-1, N, N)
        Log transition matrices.
    log_start : ndarray (N,)
        Log start probabilities.

    Returns
    -------
    states : ndarray (T,) of int
        Most likely state sequence.
    """
    T, N = log_emit.shape
    delta = np.empty((T, N), dtype=np.float64)
    psi = np.empty((T, N), dtype=np.intp)

    delta[0] = log_start + log_emit[0]

    for t in range(1, T):
        # For each state j, find best predecessor i
        # scores[i, j] = delta[t-1, i] + log_trans[t-1, i, j]
        scores = delta[t - 1, :, np.newaxis] + log_trans[t - 1]  # (N, N)
        psi[t] = scores.argmax(axis=0)
        delta[t] = log_emit[t] + scores.max(axis=0)

    # Backtrack
    states = np.empty(T, dtype=np.intp)
    states[-1] = delta[-1].argmax()
    for t in range(T - 2, -1, -1):
        states[t] = psi[t + 1, states[t + 1]]

    return states


# --- Observation preparation ---


def prepare_observations(
    cnarr: CNA,
    variants: VariantArray | None,
) -> tuple[
    list[NDArray[np.float64]],
    list[NDArray[np.float64] | None],
    list[NDArray[np.float64] | None],
    list[NDArray[np.float64]],
    list[str],
]:
    """Prepare per-arm observations for HMM.

    Returns
    -------
    log2_list : list of ndarray
        Log2 ratios per arm.
    minor_count_list : list of (ndarray or None)
        Aggregated minor allele counts per bin per arm, or None if no variants.
    depth_list : list of (ndarray or None)
        Aggregated read depths per bin per arm, or None if no variants.
    dist_list : list of ndarray
        Inter-bin distances per arm (length T-1 each).
    arm_labels : list of str
        Arm identifier strings.
    """
    log2_list: list[NDArray[np.float64]] = []
    minor_count_list: list[NDArray[np.float64] | None] = []
    depth_list: list[NDArray[np.float64] | None] = []
    dist_list: list[NDArray[np.float64]] = []
    arm_labels: list[str] = []
    _warned_no_counts = False

    for arm_label, arm in cnarr.by_arm():
        if len(arm) == 0:
            continue

        log2_obs = arm["log2"].to_numpy(dtype=np.float64)
        log2_list.append(log2_obs)
        arm_labels.append(str(arm_label))

        # Inter-bin distances
        starts = arm["start"].to_numpy(dtype=np.float64)
        ends = arm["end"].to_numpy(dtype=np.float64)
        if len(starts) > 1:
            distances = np.maximum(0.0, starts[1:] - ends[:-1])
        else:
            distances = np.array([], dtype=np.float64)
        dist_list.append(distances)

        # Allele counts per bin
        if variants is not None:
            counts = variants.baf_counts_by_ranges(arm)
            if counts is not None:
                mc = counts[0].to_numpy(dtype=np.float64)
                dp = counts[1].to_numpy(dtype=np.float64)
                minor_count_list.append(mc)
                depth_list.append(dp)
            else:
                if not _warned_no_counts:
                    logging.warning(
                        "VCF lacks alt_count/depth columns; "
                        "BAF emission disabled (log2 only)"
                    )
                    _warned_no_counts = True
                minor_count_list.append(None)
                depth_list.append(None)
        else:
            minor_count_list.append(None)
            depth_list.append(None)

    return log2_list, minor_count_list, depth_list, dist_list, arm_labels


# --- Start probabilities ---


def compute_start_probs(
    states: list[tuple[int, int | None]],
) -> NDArray[np.float64]:
    """Compute log start probabilities favoring CN=2.

    Gaussian-like weighting centered on CN=2.
    """
    cn_arr = np.array([s[0] for s in states], dtype=np.float64)
    weights = np.exp(-0.5 * ((cn_arr - 2.0) / 1.0) ** 2)
    weights /= weights.sum()
    return np.log(np.maximum(weights, 1e-300))  # type: ignore[no-any-return]


# --- Grid search ---


def _batched_emission_log2(
    log2_obs: NDArray[np.float64],
    cn_arr: NDArray[np.float64],
    purities: NDArray[np.float64],
    ploidies: NDArray[np.float64],
    log2_stdev: float,
) -> NDArray[np.float64]:
    """Compute log2 emissions for all grid points at once.

    Parameters
    ----------
    log2_obs : ndarray (T,)
    cn_arr : ndarray (N,)
    purities : ndarray (G_p,)
    ploidies : ndarray (G_l,)
    log2_stdev : float

    Returns
    -------
    ndarray (G, T, N) where G = G_p * G_l
    """
    # Expected log2 for all (purity, ploidy, cn) combinations
    # purity: (G_p, 1), ploidy: (G_l, 1), cn: (N,)
    p = purities[:, np.newaxis, np.newaxis]  # (G_p, 1, 1)
    plo = ploidies[np.newaxis, :, np.newaxis]  # (1, G_l, 1)
    cn = cn_arr[np.newaxis, np.newaxis, :]  # (1, 1, N)

    numerator = p * cn + (1.0 - p) * 2.0
    denominator = p * plo + (1.0 - p) * 2.0
    ratio = np.maximum(numerator, 1e-4) / denominator
    expected_l2 = np.log2(ratio)  # (G_p, G_l, N)

    # Reshape to (G, N) where G = G_p * G_l
    G_p, G_l, N = expected_l2.shape
    expected_l2 = expected_l2.reshape(G_p * G_l, N)  # (G, N)

    # Compute log emission: N(obs | expected, stdev)
    # log2_obs: (T,), expected_l2: (G, N) -> diff: (G, T, N)
    diff = log2_obs[np.newaxis, :, np.newaxis] - expected_l2[:, np.newaxis, :]
    log_p = -0.5 * (diff / log2_stdev) ** 2 - np.log(log2_stdev * np.sqrt(2.0 * np.pi))
    return log_p  # type: ignore[no-any-return]  # (G, T, N)


def _batched_emission_baf(
    minor_counts: NDArray[np.float64],
    depths: NDArray[np.float64],
    cn_arr: NDArray[np.float64],
    minor_arr: NDArray[np.float64],
    purities: NDArray[np.float64],
    n_ploidies: int,
    rho: float,
) -> NDArray[np.float64]:
    """Compute beta-binomial BAF emissions for all grid points at once.

    Returns
    -------
    ndarray (G, T, N)
    """
    p = purities[:, np.newaxis]  # (G_p, 1)
    cn = cn_arr[np.newaxis, :]  # (1, N)
    mn = minor_arr[np.newaxis, :]  # (1, N)

    numerator = p * mn + (1.0 - p) * 1.0
    denominator = p * cn + (1.0 - p) * 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        expected_freq = np.where(denominator > 0, numerator / denominator, 0.5)
    # (G_p, N)

    # Repeat for each ploidy (BAF doesn't depend on ploidy)
    expected_freq = np.repeat(expected_freq, n_ploidies, axis=0)  # (G, N)

    # Beta-binomial parameters with floor
    alpha = np.maximum(expected_freq * rho, 0.5)  # (G, N)
    beta = np.maximum((1.0 - expected_freq) * rho, 0.5)  # (G, N)

    k = minor_counts[np.newaxis, :, np.newaxis]  # (1, T, 1)
    n = depths[np.newaxis, :, np.newaxis]  # (1, T, 1)
    a = alpha[:, np.newaxis, :]  # (G, 1, N)
    b = beta[:, np.newaxis, :]  # (G, 1, N)

    log_baf = _betabinom_logpmf(k, n, a, b)  # (G, T, N)
    # Zero-depth bins contribute nothing
    zero_depth = depths == 0
    log_baf[:, zero_depth, :] = 0.0
    return log_baf  # type: ignore[no-any-return]  # (G, T, N)


def _forward_ll_batched(
    log_emit: NDArray[np.float64],
    log_trans: NDArray[np.float64],
    log_start: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Batched forward algorithm: compute log-likelihood for G grid points.

    Parameters
    ----------
    log_emit : ndarray (G, T, N)
        Log emission probabilities for each grid point.
    log_trans : ndarray (T-1, N, N)
        Log transition matrices (shared across grid points).
    log_start : ndarray (N,)
        Log start probabilities.

    Returns
    -------
    ndarray (G,)
        Log-likelihood for each grid point.
    """
    G, T, N = log_emit.shape
    # alpha: (G, N) — only keep current row
    alpha = log_start[np.newaxis, :] + log_emit[:, 0, :]  # (G, N)

    for t in range(1, T):
        # alpha[:, :, np.newaxis]: (G, N, 1)
        # log_trans[t-1]: (N, N) -> broadcast to (1, N, N)
        # scores: (G, N, N)
        scores = alpha[:, :, np.newaxis] + log_trans[t - 1][np.newaxis, :, :]
        # logsumexp over axis=1 (the "from" state): (G, N)
        alpha = log_emit[:, t, :] + logsumexp(scores, axis=1)

    # Final log-likelihood: logsumexp over states
    return logsumexp(alpha, axis=1)  # type: ignore[no-any-return]  # (G,)


def grid_search_purity_ploidy(
    log2_list: list[NDArray[np.float64]],
    minor_count_list: list[NDArray[np.float64] | None],
    depth_list: list[NDArray[np.float64] | None],
    log_trans_list: list[NDArray[np.float64]],
    log2_stdev: float,
    states: list[tuple[int, int | None]],
    purity_range: tuple[float, float, float],
    ploidy_range: tuple[float, float, float],
    log_start: NDArray[np.float64],
) -> pd.DataFrame:
    """Grid search over (purity, ploidy) using batched forward log-likelihood.

    All grid points are processed simultaneously per arm, with emissions
    batched across the grid and transitions shared (since they are
    distance-dependent but purity/ploidy-independent).

    Parameters
    ----------
    log2_list : list of ndarray
        Log2 observations per arm (autosomal only).
    minor_count_list : list of (ndarray or None)
        Minor allele counts per arm.
    depth_list : list of (ndarray or None)
        Read depths per arm.
    log_trans_list : list of ndarray
        Precomputed log transition matrices per arm.
    log2_stdev : float
        Emission stdev for log2.
    states : list of (cn, minor|None)
        State definitions.
    purity_range : (min, max, step)
        Purity grid parameters.
    ploidy_range : (min, max, step)
        Ploidy grid parameters.
    log_start : ndarray (N,)
        Log start probabilities.

    Returns
    -------
    DataFrame
        Columns: purity, ploidy, score. Sorted by score descending.
    """
    include_baf = any(mc is not None for mc in minor_count_list)

    purities = np.arange(purity_range[0], purity_range[1] + 1e-9, purity_range[2])
    ploidies = np.arange(ploidy_range[0], ploidy_range[1] + 1e-9, ploidy_range[2])
    G = len(purities) * len(ploidies)

    cn_arr = np.array([s[0] for s in states], dtype=np.float64)
    minor_arr = (
        np.array([s[1] if s[1] is not None else 0 for s in states], dtype=np.float64)
        if include_baf
        else None
    )

    total_ll = np.zeros(G, dtype=np.float64)

    for arm_idx in range(len(log2_list)):
        l2 = log2_list[arm_idx]
        lt = log_trans_list[arm_idx]
        if len(l2) == 0:
            continue

        # Compute batched emissions: (G, T, N)
        log_emit = _batched_emission_log2(l2, cn_arr, purities, ploidies, log2_stdev)

        if (
            include_baf
            and minor_count_list[arm_idx] is not None
            and depth_list[arm_idx] is not None
        ):
            assert minor_arr is not None
            log_emit += _batched_emission_baf(
                minor_count_list[arm_idx],  # type: ignore[arg-type]
                depth_list[arm_idx],  # type: ignore[arg-type]
                cn_arr,
                minor_arr,
                purities,
                len(ploidies),
                rho=BETABINOM_RHO,
            )

        if len(l2) == 1:
            total_ll += logsumexp(log_start[np.newaxis, :] + log_emit[:, 0, :], axis=1)
            continue

        total_ll += _forward_ll_batched(log_emit, lt, log_start)

    # Build result grid
    pur_grid = np.repeat(purities, len(ploidies))
    plo_grid = np.tile(ploidies, len(purities))
    df = pd.DataFrame({"purity": pur_grid, "ploidy": plo_grid, "score": total_ll})
    return df.sort_values("score", ascending=False).reset_index(drop=True)


# --- Method parameters ---


def _method_params(method: str) -> dict:
    """Get HMM parameters for the given method variant."""
    match method:
        case "hmm-germline":
            return {"max_copies": 4, "fix_purity": 1.0, "fix_ploidy": 2.0}
        case "hmm-tumor":
            return {"max_copies": 8, "min_purity": 0.1}
        case "hmm":
            return {"max_copies": 6, "min_purity": 0.2}
        case _:
            raise ValueError(f"Unknown HMM method: {method!r}")


# --- Main entry point ---


def segment_hmm(
    cnarr: CNA,
    method: str,
    diploid_parx_genome: str | None,
    window: int | None = None,
    variants: VariantArray | None = None,
    processes: int = 1,
) -> CNA:
    """Segment bins using Hidden Markov Model with Viterbi decoding.

    Uses a pure numpy/scipy HMM with distance-dependent transitions,
    Gaussian log2 ratio emissions, and optional beta-binomial BAF emissions
    on aggregated allele counts — parameterized by (purity, ploidy).
    For somatic methods ('hmm', 'hmm-tumor'), purity and ploidy are estimated
    via grid search over marginal likelihood on autosomal arms.

    Parameters
    ----------
    cnarr : CopyNumArray
        Bin-level copy ratios (.cnr file) to segment.
    method : str
        HMM variant: 'hmm', 'hmm-tumor', or 'hmm-germline'.
    diploid_parx_genome : str or None
        Reference genome name for pseudo-autosomal region handling.
    window : int or float or None
        Smoothing window width.
    variants : VariantArray or None
        Variant allele frequency data (from VCF) for joint segmentation.
    processes : int
        Number of parallel processes (reserved for future use).

    Returns
    -------
    CopyNumArray
        Segmented copy number data (.cns format).
    """
    params = _method_params(method)
    include_baf = variants is not None
    max_copies: int = params["max_copies"]

    # 1. Smooth log2
    orig_log2 = cnarr["log2"].to_numpy().copy()
    cnarr["log2"] = cnarr.smooth_log2(window)

    # 2. Estimate stdev robustly
    all_log2 = cnarr["log2"].to_numpy()
    log2_stdev = max(float(biweight_midvariance(all_log2, initial=0)), LOG2_STDEV_FLOOR)

    # 3. Enumerate states
    states = enumerate_states(max_copies, include_baf)
    n_states = len(states)
    cn_values = np.array([s[0] for s in states], dtype=np.intp)
    log_start = compute_start_probs(states)

    # 4. Prepare observations per arm
    log2_list, minor_count_list, depth_list, dist_list, arm_labels = (
        prepare_observations(cnarr, variants)
    )
    logging.info(
        "HMM segmentation: %d arms, %d states, method=%s",
        len(log2_list),
        n_states,
        method,
    )

    # 5. Precompute transition matrices per arm (distance-dependent, shared across grid)
    log_trans_list: list[NDArray[np.float64]] = []
    for dists in dist_list:
        if len(dists) > 0:
            lt = log_transition_matrices(dists, n_states, cn_values)
        else:
            lt = np.empty((0, n_states, n_states), dtype=np.float64)
        log_trans_list.append(lt)

    # 6. Choose purity/ploidy
    purity_results: pd.DataFrame | None = None
    if "fix_purity" in params:
        # hmm-germline: fixed purity and ploidy, no grid search
        best_purity = params["fix_purity"]
        best_ploidy = params["fix_ploidy"]
        logging.info("Using fixed purity=%.2f, ploidy=%.1f", best_purity, best_ploidy)
    else:
        # Grid search on autosomal arms
        auto_cnarr = cnarr.autosomes(diploid_parx_genome=diploid_parx_genome)
        auto_arms = {str(label) for label, _ in auto_cnarr.by_arm()}

        auto_idx = [i for i, lab in enumerate(arm_labels) if lab in auto_arms]
        auto_log2 = [log2_list[i] for i in auto_idx]
        auto_mc = [minor_count_list[i] for i in auto_idx]
        auto_dp = [depth_list[i] for i in auto_idx]
        auto_trans = [log_trans_list[i] for i in auto_idx]

        if not auto_log2:
            logging.warning(
                "No autosomal chromosomes found; "
                "skipping purity/ploidy estimation, using defaults"
            )
            best_purity = 1.0
            best_ploidy = 2.0
        else:
            min_pur = params.get("min_purity", 0.2)
            purity_results = grid_search_purity_ploidy(
                auto_log2,
                auto_mc,
                auto_dp,
                auto_trans,
                log2_stdev,
                states,
                purity_range=(min_pur, 1.0, 0.05),
                ploidy_range=(1.5, 5.0, 0.5),
                log_start=log_start,
            )
            best_purity = float(purity_results.iloc[0]["purity"])
            best_ploidy = float(purity_results.iloc[0]["ploidy"])
            logging.info(
                "Best purity=%.2f, ploidy=%.1f (score=%.1f)",
                best_purity,
                best_ploidy,
                float(purity_results.iloc[0]["score"]),
            )

    # 7. Run Viterbi on all arms with best purity/ploidy
    all_states: list[NDArray[np.int_]] = []
    for arm_idx in range(len(log2_list)):
        l2 = log2_list[arm_idx]
        mc = minor_count_list[arm_idx]
        dp = depth_list[arm_idx]
        lt = log_trans_list[arm_idx]

        if len(l2) == 0:
            all_states.append(np.array([], dtype=np.intp))
            continue

        log_emit = log_emission_probs(
            l2, mc, dp, states, best_purity, best_ploidy, log2_stdev
        )

        if len(l2) == 1:
            # Single-bin arm: pick best state directly
            best = int((log_start + log_emit[0]).argmax())
            all_states.append(np.array([best], dtype=np.intp))
            continue

        arm_states = viterbi_log(log_emit, lt, log_start)
        all_states.append(arm_states)

    # 8. Map state indices to CN values
    state_indices = np.concatenate(all_states)

    logging.info("Predicted %d state values", len(state_indices))
    logging.debug(
        "State distribution: %s", np.bincount(state_indices, minlength=n_states)
    )

    # 9. Restore original log2, squash
    cnarr["log2"] = orig_log2
    cnarr["probes"] = 1
    segarr = squash_by_groups(
        cnarr, pd.Series(state_indices, index=cnarr.data.index), by_arm=True
    )
    if not (segarr.start < segarr.end).all():
        bad_segs = segarr[segarr.start >= segarr.end]
        logging.warning("Bad segments:\n%s", bad_segs.data)

    # 10. Store grid results for optional output
    if purity_results is not None:
        segarr.meta["purity_results"] = purity_results

    return segarr


# --- BAF re-segmentation (for non-HMM methods) ---


def variants_in_segment(varr, segment, min_variants=MIN_VARIANTS_THRESHOLD):
    """Re-segment a segment based on variant allele frequencies.

    Uses a simple 2-state HMM (neutral BAF=0.5 vs alt BAF=0.67) with
    homogeneous transitions and Viterbi decoding.

    Called for non-HMM segmentation methods to refine breakpoints using BAF.
    """
    results = None

    if len(varr) > min_variants:
        observations = varr.mirrored_baf(above_half=True)

        if len(observations) > 0:
            obs = observations.to_numpy(dtype=np.float64)

            # 2-state model: neutral (BAF~0.5) and alt (BAF~0.67)
            means = np.array([0.5, 0.67])
            stdev = 0.1

            # Log emission: N(obs | mean, stdev)
            diff = obs[:, np.newaxis] - means[np.newaxis, :]  # (T, 2)
            log_emit = -0.5 * (diff / stdev) ** 2 - np.log(stdev * np.sqrt(2.0 * np.pi))

            # Homogeneous transitions: strongly prefer staying
            p_stay = 0.99
            trans_2x2 = np.array([[p_stay, 1.0 - p_stay], [1.0 - p_stay, p_stay]])
            log_trans_2x2 = np.log(trans_2x2)
            # Broadcast to (T-1, 2, 2)
            T = len(obs)
            log_trans = np.broadcast_to(log_trans_2x2, (max(T - 1, 0), 2, 2)).copy()

            # Start: prefer neutral
            log_start = np.log(np.array([0.95, 0.05]))

            if T > 1:
                state_seq = viterbi_log(log_emit, log_trans, log_start)
            else:
                state_seq = np.array(
                    [int((log_start + log_emit[0]).argmax())], dtype=np.intp
                )

            # Merge adjacent bins with the same state to create segments
            fake_cnarr = CNA(varr.add_columns(weight=1, log2=0, gene=".").data)
            results = squash_by_groups(
                fake_cnarr, varr.as_series(state_seq), by_arm=False
            )
            assert (results.start < results.end).all()

    if results is not None and len(results) > 1:
        logging.info(
            "Segment %s:%d-%d on allele freqs for %d additional breakpoints",
            segment.chromosome,
            segment.start,
            segment.end,
            len(results) - 1,
        )
        # Place breakpoints midway between SNVs
        mid_breakpoints = (
            results.start.to_numpy()[1:] + results.end.to_numpy()[:-1]
        ) // 2
        starts = np.concatenate([[segment.start], mid_breakpoints])
        ends = np.concatenate([mid_breakpoints, [segment.end]])
        dframe = pd.DataFrame(
            {
                "chromosome": segment.chromosome,
                "start": starts,
                "end": ends,
                "gene": segment.gene,
                "log2": segment.log2,
                "probes": results["probes"],
            }
        )
        bad_segs_idx = dframe.start >= dframe.end
        if bad_segs_idx.any():
            raise RuntimeError(
                f"Improper post-processing of segment {segment} -- "
                f"{bad_segs_idx.sum()} bins start >= end:\n{dframe[bad_segs_idx]}\n"
            )

    else:
        dframe = pd.DataFrame(
            {
                "chromosome": segment.chromosome,
                "start": segment.start,
                "end": segment.end,
                "gene": segment.gene,
                "log2": segment.log2,
                "probes": segment.probes,
            },
            index=[0],
        )

    return dframe
