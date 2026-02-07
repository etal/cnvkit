"""Segmentation by Hidden Markov Model."""

from __future__ import annotations

import collections
import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import scipy.special

from ..cnary import CopyNumArray as CNA
from ..descriptives import biweight_midvariance
from ..segfilters import squash_by_groups

if TYPE_CHECKING:
    from cnvlib.vary import VariantArray
    from numpy import ndarray
    from pomegranate.hmm.dense_hmm import DenseHMM

# Configuration constants
DEFAULT_MAX_ITER = 100  # Reduced from 100000 for performance with PyTorch backend
DEFAULT_TOLERANCE = 1e-4
DEFAULT_INERTIA = 0.1
TRANSITION_SELF_PREFERENCE = 100  # Weight for staying in same state
MIN_VARIANTS_THRESHOLD = 50

# Optional pomegranate import for HMM functionality
try:
    import pomegranate
    import pomegranate.distributions
    from pomegranate.distributions import Normal
    from pomegranate.hmm import DenseHMM

    # Version compatibility check
    POMEGRANATE_VERSION = tuple(map(int, pomegranate.__version__.split(".")[:2]))
    if POMEGRANATE_VERSION < (1, 0):
        raise ImportError(
            f"pomegranate >= 1.0.0 required, found {pomegranate.__version__}"
        )

    HMM_AVAILABLE = True
except ImportError as e:
    HMM_AVAILABLE = False
    HMM_IMPORT_ERROR = str(e)


def segment_hmm(
    cnarr: CNA,
    method: str,
    diploid_parx_genome: str | None,
    window: Any | None = None,
    variants: VariantArray | None = None,
    processes: int = 1,
) -> CNA:
    """Segment bins using Hidden Markov Model (HMM) with Viterbi algorithm.

    Applies probabilistic HMM-based segmentation to infer discrete copy number
    states from continuous log2 ratio data. The Viterbi algorithm finds the
    most likely sequence of hidden states (copy number levels) given the
    observed coverage data.

    This method can model 3-state (loss/neutral/gain) or 5-state (including
    deep deletions and amplifications) copy number profiles, with either
    flexible or fixed state means.

    Parameters
    ----------
    cnarr : CopyNumArray
        Bin-level copy ratios (.cnr file) to segment. Should be normalized
        and centered around log2=0 for neutral copy number.
    method : str
        HMM variant to use. Options:

        - 'hmm': 3-state model (loss/neutral/gain) with flexible means
        - 'hmm-tumor': 5-state model (deep loss/loss/neutral/gain/amplification)
          with flexible means, optimized for tumor samples
        - 'hmm-germline': 3-state model with fixed means, optimized for
          germline samples
    diploid_parx_genome : str, optional
        Reference genome name (e.g., 'hg19', 'hg38') for pseudo-autosomal
        region handling. Affects modeling of X/Y chromosomes.
        Default: None
    window : int, optional
        Smoothing window width (number of bins) applied before HMM fitting.
        If None, uses weighted smoothing based on bin statistics.
        Default: None
    variants : VariantArray, optional
        Variant allele frequency data (from VCF). Currently not actively used
        in HMM segmentation but reserved for future joint segmentation of
        copy ratios and allele frequencies.
        Default: None
    processes : int, optional
        Number of parallel processes for model building (when processing
        multiple chromosomes). Currently limited to 1 due to pomegranate
        constraints.
        Default: 1

    Returns
    -------
    CopyNumArray
        Segmented copy number data (.cns file) with discrete states.
        Each segment represents a region with homogeneous copy number state.
        Contains columns: chromosome, start, end, gene, log2, probes.

    Raises
    ------
    ImportError
        If pomegranate library (>=1.0.0) is not installed. Install with:
        ``pip install pomegranate>=1.0.0``
    ValueError
        If no valid observations are available for HMM prediction (e.g.,
        all-zero coverage).

    Notes
    -----
    The HMM segmentation process:

    1. **Smoothing**: Apply weighted smoothing to log2 ratios to reduce noise
       while preserving biological signal.

    2. **Model building**: Construct HMM with discrete states representing
       copy number levels. State means and transition probabilities are
       learned from the data (except for 'hmm-germline').

    3. **Viterbi decoding**: Find the most likely sequence of states that
       explains the observed log2 ratios.

    4. **Segmentation**: Merge adjacent bins with the same predicted state
       into contiguous segments.

    **State models**:

    - 3-state (hmm, hmm-germline): deletion (CN<2), neutral (CN=2), gain (CN>2)
    - 5-state (hmm-tumor): deep deletion (CN=0), loss (CN=1), neutral (CN=2),
      gain (CN=3), amplification (CN≥4)

    **Advantages over CBS**:

    - Probabilistic framework with explicit copy number states
    - Can handle noisy data more robustly
    - Joint modeling of multiple features (future enhancement)

    **Limitations**:

    - Requires pomegranate library (~500MB+ with PyTorch dependencies)
    - Slower than CBS for large datasets
    - Less established in CNV literature than CBS

    See Also
    --------
    hmm_get_model : Builds the HMM from observations
    as_observation_matrix : Converts copy ratios to HMM observation format
    segment_cbs : Alternative segmentation using Circular Binary Segmentation

    Examples
    --------
    Basic 3-state HMM segmentation:
    >>> segments = segment_hmm(
    ...     cnarr=cnr,
    ...     method='hmm',
    ...     diploid_parx_genome='hg38'
    ... )

    5-state tumor-specific segmentation:
    >>> segments = segment_hmm(
    ...     cnarr=tumor_cnr,
    ...     method='hmm-tumor',
    ...     diploid_parx_genome='hg38',
    ...     window=50  # Custom smoothing window
    ... )

    Germline segmentation with fixed state means:
    >>> segments = segment_hmm(
    ...     cnarr=germline_cnr,
    ...     method='hmm-germline',
    ...     diploid_parx_genome=None
    ... )
    """
    if not HMM_AVAILABLE:
        raise ImportError(
            f"HMM segmentation requires pomegranate >= 1.0.0. "
            f"Install with: pip install pomegranate>=1.0.0. Error: {HMM_IMPORT_ERROR}"
        )
    # NB: Incorporate weights into smoothed log2 estimates
    # (Useful kludge until weighted HMM is in place)
    orig_log2 = cnarr["log2"].to_numpy().copy()
    cnarr["log2"] = cnarr.smooth_log2(window)

    logging.info("Building model from observations")
    model = hmm_get_model(cnarr, method, diploid_parx_genome, processes)

    logging.info("Predicting states from model")
    observations = as_observation_matrix(cnarr)
    # Format observations for prediction (same as in fit)
    # Filter out empty observations that would cause tensor shape issues
    valid_observations = [obs for obs in observations if len(obs) > 0]

    if not valid_observations:
        raise ValueError("No valid observations available for HMM prediction")

    # Convert to 3D tensor format expected by pomegranate 1.0: (n_sequences, max_length, n_features)
    max_length = max(len(obs) for obs in valid_observations)
    n_sequences = len(valid_observations)
    n_features = 1

    # Create padded 3D tensor with explicit float32 dtype for PyTorch compatibility
    X_3d = np.zeros((n_sequences, max_length, n_features), dtype=np.float32)
    sequence_lengths = []

    for i, obs in enumerate(valid_observations):
        seq_len = len(obs)
        X_3d[i, :seq_len, 0] = obs
        sequence_lengths.append(seq_len)

    # Get predictions for all sequences at once
    all_predictions = model.predict(X_3d)

    # Extract only the valid parts of each sequence (removing padding)
    states = np.concatenate(
        [all_predictions[i, : sequence_lengths[i]] for i in range(n_sequences)]
    )

    logging.info("Done, now finalizing")
    logging.debug("Predicted states: %s", states[:100])
    logging.debug(str(collections.Counter(states)))
    logging.debug("Observations: %s", observations[0][:100])
    logging.debug("Edges: %s", model.edges)

    # Merge adjacent bins with the same state to create segments
    cnarr["log2"] = orig_log2
    cnarr["probes"] = 1
    segarr = squash_by_groups(
        cnarr, pd.Series(states, index=cnarr.data.index), by_arm=True
    )
    if not (segarr.start < segarr.end).all():
        bad_segs = segarr[segarr.start >= segarr.end]
        logging.warning("Bad segments:\n%s", bad_segs.data)
    return segarr


def hmm_get_model(
    cnarr: CNA, method: str, diploid_parx_genome: str | None, processes: int
) -> pomegranate.hmm.dense_hmm.DenseHMM:
    """

    Parameters
    ----------
    cnarr : CopyNumArray
        The normalized bin-level values to be segmented.
    method : string
        One of 'hmm', 'hmm-tumor', 'hmm-germline'.
    diploid_parx_genome : string
        Whether to include PAR1/2 from chr X within the autosomes.
    processes : int
        Number of parallel jobs to run.

    Returns
    -------
    model :
        A pomegranate HiddenMarkovModel trained on the given dataset.
    """
    assert method in ("hmm-tumor", "hmm-germline", "hmm")
    observations = as_observation_matrix(
        cnarr.autosomes(diploid_parx_genome=diploid_parx_genome)
    )

    # Estimate standard deviation from the full distribution, robustly
    stdev = biweight_midvariance(np.concatenate(observations), initial=0)
    if method == "hmm-germline":
        state_names = ["loss", "neutral", "gain"]
        distributions = [
            Normal([-1.0], [stdev], covariance_type="diag", frozen=True, inertia=0.8),
            Normal([0.0], [stdev], covariance_type="diag", frozen=True, inertia=0.8),
            Normal([0.585], [stdev], covariance_type="diag", frozen=True, inertia=0.8),
        ]
    elif method == "hmm-tumor":
        state_names = ["del", "loss", "neutral", "gain", "amp"]
        distributions = [
            Normal([-2.0], [stdev], covariance_type="diag", frozen=False, inertia=0.8),
            Normal([-0.5], [stdev], covariance_type="diag", frozen=False, inertia=0.8),
            Normal([0.0], [stdev], covariance_type="diag", frozen=True, inertia=0.8),
            Normal([0.3], [stdev], covariance_type="diag", frozen=False, inertia=0.8),
            Normal([1.0], [stdev], covariance_type="diag", frozen=False, inertia=0.8),
        ]
    else:
        state_names = ["loss", "neutral", "gain"]
        distributions = [
            Normal([-1.0], [stdev], covariance_type="diag", frozen=False, inertia=0.8),
            Normal([0.0], [stdev], covariance_type="diag", frozen=False, inertia=0.8),
            Normal([0.585], [stdev], covariance_type="diag", frozen=False, inertia=0.8),
        ]

    n_states = len(distributions)
    # Starts -- prefer neutral
    binom_coefs = scipy.special.binom(n_states - 1, range(n_states))
    start_probabilities = (binom_coefs / binom_coefs.sum()).astype(np.float32)

    # Prefer to keep the current state in each transition
    # All other transitions are equally likely, to start
    transition_matrix = (
        np.identity(n_states) * TRANSITION_SELF_PREFERENCE
        + np.ones((n_states, n_states)) / n_states
    )
    # Rescale so max is 1.0
    transition_matrix = (transition_matrix / transition_matrix.max()).astype(np.float32)

    model = DenseHMM(
        distributions=distributions,
        edges=transition_matrix,
        starts=start_probabilities,
        ends=start_probabilities,
        inertia=DEFAULT_INERTIA,  # edge_inertia from old API
        verbose=False,
        max_iter=DEFAULT_MAX_ITER,
        tol=DEFAULT_TOLERANCE,
    )

    # Convert observations to proper multidimensional format for pomegranate 1.0
    # Filter out empty observations that would cause tensor shape issues
    valid_observations = [obs for obs in observations if len(obs) > 0]

    if not valid_observations:
        raise ValueError("No valid observations available for HMM training")

    # Convert to 3D tensor format expected by pomegranate 1.0: (n_sequences, max_length, n_features)
    max_length = max(len(obs) for obs in valid_observations)
    n_sequences = len(valid_observations)
    n_features = 1

    # Create padded 3D tensor with explicit float32 dtype for PyTorch compatibility
    X_3d = np.zeros((n_sequences, max_length, n_features), dtype=np.float32)
    sequence_lengths = []

    for i, obs in enumerate(valid_observations):
        seq_len = len(obs)
        X_3d[i, :seq_len, 0] = obs
        sequence_lengths.append(seq_len)

    model.fit(X=X_3d, sample_weight=sequence_lengths)
    return model


def as_observation_matrix(cnarr: CNA, variants: Any | None = None) -> list[ndarray]:
    """Extract HMM fitting values from `cnarr`.

    For each chromosome arm, extract log2 ratios as a numpy array.

    Future: If VCF of variants is given, or 'baf' column has already been
    added to `cnarr` from the same, then the BAF values are a second row/column
    in each numpy array.

    Returns: List of numpy.ndarray, one per chromosome arm.
    """
    # TODO incorporate weights -- currently handled by smoothing
    # TODO incorporate inter-bin distances
    observations = [arm.log2.to_numpy() for _c, arm in cnarr.by_arm()]
    return observations


def variants_in_segment(varr, segment, min_variants=MIN_VARIANTS_THRESHOLD):
    results = None

    if len(varr) > min_variants:
        observations = varr.mirrored_baf(above_half=True)

        # Check if observations are valid (not empty) before proceeding with HMM
        if len(observations) > 0:
            state_names = ["neutral", "alt"]
            distributions = [
                Normal([0.5], [0.1], covariance_type="diag", frozen=True, inertia=0.8),
                Normal([0.67], [0.1], covariance_type="diag", frozen=True, inertia=0.8),
            ]
            n_states = len(distributions)
            # Starts -- prefer neutral
            start_probabilities = np.array([0.95, 0.05], dtype=np.float32)
            # Prefer to keep the current state in each transition
            # All other transitions are equally likely, to start
            transition_matrix = (
                np.identity(n_states) * TRANSITION_SELF_PREFERENCE
                + np.ones((n_states, n_states)) / n_states
            )
            # Rescale so max is 1.0
            transition_matrix = (transition_matrix / transition_matrix.max()).astype(
                np.float32
            )
            model = DenseHMM(
                distributions=distributions,
                edges=transition_matrix,
                starts=start_probabilities,
                ends=start_probabilities,
                inertia=DEFAULT_INERTIA,  # edge_inertia from old API
                verbose=False,
                max_iter=DEFAULT_MAX_ITER,
                tol=DEFAULT_TOLERANCE,
            )

            # Convert observations to 3D tensor format for pomegranate 1.0
            # Shape: (n_sequences=1, sequence_length, n_features=1)
            X_3d = observations.reshape(1, -1, 1).astype(np.float32)

            model.fit(X=X_3d, sample_weight=[len(observations)])
            states = model.predict(X_3d)[0]  # Extract the single sequence

            logging.info("Done, now finalizing")
            logging.debug("Predicted states: %s", states[:100])
            logging.debug(str(collections.Counter(states)))
            logging.debug("Edges: %s", model.edges)

            # Merge adjacent bins with the same state to create segments
            fake_cnarr = CNA(varr.add_columns(weight=1, log2=0, gene=".").data)
            results = squash_by_groups(fake_cnarr, varr.as_series(states), by_arm=False)
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
        # XXX TODO use original cnarr bin boundaries to select/adjust breakpoint
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
                # 'baf': results['mean'],
                "gene": segment.gene,  # '-'
                "log2": segment.log2,
                "probes": results["probes"],
                # 'weight': (segment.weight * results['probes']
                #            / (segment.end - segment.start)),
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
                "gene": segment.gene,  # '-',
                "log2": segment.log2,
                "probes": segment.probes,
                # 'weight': segment.weight,
            },
            index=[0],
        )

    return dframe
