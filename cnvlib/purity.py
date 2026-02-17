"""Estimate tumor purity and ploidy from segment copy ratios.

Uses a Gaussian KDE grid search over (purity, ploidy) parameter space,
scoring each candidate by evaluating the KDE at expected copy-number peak
positions. Optionally integrates BAF from a VCF to sharpen the estimate.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

if TYPE_CHECKING:
    from cnvlib.cnary import CopyNumArray
    from cnvlib.vary import VariantArray


def do_purity(
    segments: CopyNumArray,
    variants: VariantArray | None = None,
    min_purity: float = 0.1,
    max_purity: float = 1.0,
    purity_step: float = 0.01,
    min_ploidy: float = 1.5,
    max_ploidy: float = 5.0,
    ploidy_step: float = 0.1,
    max_copies: int = 8,
    baf_weight: float = 1.0,
    method: str = "kde",
) -> pd.DataFrame:
    """Estimate tumor purity and ploidy by grid search.

    Score each (purity, ploidy) combination using a Gaussian KDE fit to segment
    copy ratios, optionally augmented with b-allele frequencies.

    Parameters
    ----------
    segments : CopyNumArray
        Segmented copy ratios (.cns).
    variants : VariantArray or None
        Heterozygous SNP data from a VCF, for BAF-based scoring.
    min_purity, max_purity, purity_step : float
        Range and step size for purity grid.
    min_ploidy, max_ploidy, ploidy_step : float
        Range and step size for ploidy grid.
    max_copies : int
        Maximum integer copy number to consider.
    baf_weight : float
        Weight for the BAF term relative to the ratio term.
    method : str
        Estimation method (currently only 'kde').

    Returns
    -------
    pd.DataFrame
        Columns: purity, ploidy, score; sorted by score descending.
    """
    if method != "kde":
        raise ValueError(f"Unsupported method: {method!r}")

    # Validate grid parameters
    if min_purity >= max_purity:
        raise ValueError(
            f"min_purity ({min_purity}) must be less than max_purity ({max_purity})"
        )
    if min_ploidy >= max_ploidy:
        raise ValueError(
            f"min_ploidy ({min_ploidy}) must be less than max_ploidy ({max_ploidy})"
        )
    if purity_step <= 0 or ploidy_step <= 0:
        raise ValueError("Step sizes must be positive")

    # Extract segment log2 ratios and convert to linear space
    log2_ratios = segments["log2"].to_numpy()
    linear_ratios = 2.0**log2_ratios

    # Segment weights: prefer 'weight' column, fall back to 'probes'
    if "weight" in segments:
        weights = segments["weight"].to_numpy().astype(float)
    elif "probes" in segments:
        weights = segments["probes"].to_numpy().astype(float)
    else:
        weights = np.ones(len(segments))

    # Require positive weights; drop zero-weight segments
    mask = weights > 0
    if not mask.all():
        linear_ratios = linear_ratios[mask]
        weights = weights[mask]

    if len(linear_ratios) < 3:
        logging.warning("Too few segments (%d) for KDE estimation", len(linear_ratios))
        return pd.DataFrame(columns=["purity", "ploidy", "score"])

    # Build KDE on linear ratios
    try:
        ratio_kde = gaussian_kde(linear_ratios, weights=weights)
    except np.linalg.LinAlgError:
        logging.warning(
            "KDE failed (degenerate segment ratios); cannot estimate purity"
        )
        return pd.DataFrame(columns=["purity", "ploidy", "score"])

    # Optionally build BAF KDE
    baf_kde = None
    if variants is not None:
        baf_values = variants.baf_by_ranges(segments, above_half=True)
        baf_values = baf_values.dropna().to_numpy()
        if len(baf_values) >= 3:
            try:
                baf_kde = gaussian_kde(baf_values)
            except np.linalg.LinAlgError:
                logging.warning("BAF KDE failed (degenerate data); ignoring BAF")
        else:
            logging.warning(
                "Too few BAF values (%d) for KDE; ignoring BAF", len(baf_values)
            )

    # Build grid
    purities = np.arange(min_purity, max_purity + purity_step / 2, purity_step)
    ploidies = np.arange(min_ploidy, max_ploidy + ploidy_step / 2, ploidy_step)
    logging.info(
        "Searching %d x %d purity/ploidy grid (%d points)",
        len(purities),
        len(ploidies),
        len(purities) * len(ploidies),
    )

    results = _score_grid(
        ratio_kde, baf_kde, purities, ploidies, max_copies, baf_weight
    )
    results = results.sort_values("score", ascending=False, ignore_index=True)
    return results


def _score_grid(
    ratio_kde: gaussian_kde,
    baf_kde: gaussian_kde | None,
    purities: np.ndarray,
    ploidies: np.ndarray,
    max_copies: int,
    baf_weight: float,
) -> pd.DataFrame:
    """Score each (purity, ploidy) grid point."""
    rows = []
    cn_states = np.arange(0, max_copies + 1)
    for purity in purities:
        # BAF expected values depend only on purity, not ploidy -- compute once
        baf_score = 0.0
        if baf_kde is not None:
            baf_points = []
            for total_cn in cn_states[cn_states >= 1]:
                for minor in range(total_cn // 2 + 1):
                    baf_points.append(_expected_baf(minor, total_cn, purity))
            baf_score = float(np.sum(baf_kde(baf_points)))

        for ploidy in ploidies:
            # Score from copy-ratio KDE
            expected_ratios: np.ndarray = _expected_ratio(  # type: ignore[assignment]
                cn_states, purity, ploidy
            )
            # Only evaluate at positive expected ratios
            valid: np.ndarray = expected_ratios > 0
            if not valid.any():
                rows.append((purity, ploidy, 0.0))
                continue
            ratio_score = float(np.sum(ratio_kde(expected_ratios[valid])))

            score = ratio_score + baf_weight * baf_score
            rows.append((purity, ploidy, score))
    return pd.DataFrame(rows, columns=["purity", "ploidy", "score"])


def _expected_ratio(
    cn: int | np.ndarray, purity: float, ploidy: float
) -> float | np.ndarray:
    """Expected linear copy ratio for a given CN state.

    Formula: (purity * cn + (1 - purity) * 2) / (purity * ploidy + (1 - purity) * 2)
    """
    return (purity * cn + (1 - purity) * 2) / (purity * ploidy + (1 - purity) * 2)


def _expected_baf(minor: int, total_cn: int, purity: float) -> float:
    """Expected BAF for a minor/total CN state.

    Formula: (purity * minor + (1 - purity) * 1) / (purity * total_cn + (1 - purity) * 2)
    """
    return (purity * minor + (1 - purity) * 1) / (purity * total_cn + (1 - purity) * 2)


def read_purity_tsv(fname: str) -> tuple[float, float]:
    """Read the best purity and ploidy from a purity output TSV.

    Parameters
    ----------
    fname : str
        Path to a purity output file (tab-separated with columns
        purity, ploidy, score).

    Returns
    -------
    tuple of (float, float)
        The purity and ploidy values from the first (best-scoring) row.

    Raises
    ------
    ValueError
        If purity is not in (0, 1] or ploidy is < 1.
    """
    table = pd.read_csv(fname, sep="\t", nrows=1)
    pur = float(table["purity"].iloc[0])
    plo = float(table["ploidy"].iloc[0])
    if not 0.0 < pur <= 1.0:
        raise ValueError(f"Purity in {fname} must be in (0, 1], got {pur}")
    if plo < 1:
        raise ValueError(f"Ploidy in {fname} must be >= 1, got {plo}")
    return pur, plo
