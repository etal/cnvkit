"""Robust metrics to evaluate performance of copy number estimates."""

from __future__ import annotations
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
import pandas as pd

from . import descriptives

if TYPE_CHECKING:
    from collections.abc import Iterator
    from cnvlib.cnary import CopyNumArray
    from numpy import float64, ndarray


def do_metrics(
    cnarrs: CopyNumArray,
    segments: Optional[CopyNumArray] = None,
    skip_low: bool = False,
) -> pd.DataFrame:
    """Compute coverage deviations and other metrics for self-evaluation.

    Parameters
    ----------
    cnarrs : CopyNumArray or list of CopyNumArray
        Bin-level copy number data for one or more samples.
    segments : CopyNumArray, list of CopyNumArray, or None, optional
        Segmented copy number data. If None, computes metrics without
        segment-based residuals.
    skip_low : bool, optional
        Skip bins with low coverage. Default is False.

    Returns
    -------
    pd.DataFrame
        Metrics table with columns: sample, segments, stdev, mad, iqr, bivar.
        Each row contains quality metrics for one sample.
    """
    # Catch if passed args are single CopyNumArrays instead of lists
    from .cnary import CopyNumArray as CNA

    if isinstance(cnarrs, CNA):
        cnarrs = [cnarrs]
    if isinstance(segments, CNA):
        segments = [segments]
    elif segments is None:
        segments = [None]
    else:
        segments = list(segments)
    if skip_low:
        cnarrs = (cna.drop_low_coverage() for cna in cnarrs)
    rows = (
        (
            cna.meta.get("filename", cna.sample_id),
            len(seg) if seg is not None else "-",
            *ests_of_scale(cna.residuals(seg).to_numpy()),
        )
        for cna, seg in zip_repeater(cnarrs, segments)
    )
    colnames = ["sample", "segments", "stdev", "mad", "iqr", "bivar"]
    return pd.DataFrame.from_records(rows, columns=colnames)


def zip_repeater(
    iterable: Iterator[Any], repeatable: list[CopyNumArray]
) -> Iterator[tuple[CopyNumArray, CopyNumArray]]:
    """Repeat a single segmentation to match the number of copy ratio inputs"""
    rpt_len = len(repeatable)
    if rpt_len == 1:
        rpt = repeatable[0]
        for it in iterable:
            yield it, rpt
    else:
        i = -1
        for i, (it, rpt) in enumerate(zip(iterable, repeatable, strict=False)):
            yield it, rpt
        # Require lengths to match
        if i + 1 != rpt_len:
            raise ValueError(
                "Number of unsegmented and segmented input files did not match "
                + f"({i} vs. {rpt_len})"
            )


def ests_of_scale(deviations: ndarray) -> tuple[float64, float64, float64, float64]:
    """Estimators of scale: standard deviation, MAD, biweight midvariance.

    Calculates all of these values for an array of deviations and returns them
    as a tuple.
    """
    std = np.std(deviations, dtype=np.float64)
    mad = descriptives.median_absolute_deviation(deviations)
    iqr = descriptives.interquartile_range(deviations)
    biw = descriptives.biweight_midvariance(deviations)
    return (std, mad, iqr, biw)
