"""Robust metrics to evaluate performance of copy number estimates.
"""
from __future__ import absolute_import, division, print_function
from builtins import zip

import numpy as np
import pandas as pd

from . import descriptives


def do_metrics(cnarrs, segments=None, skip_low=False):
    """Compute coverage deviations and other metrics for self-evaluation."""
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
    rows = ((cna.meta.get("filename", cna.sample_id),
             len(seg) if seg is not None else '-'
            ) + ests_of_scale(cna.residuals(seg))
            for cna, seg in zip_repeater(cnarrs, segments))
    colnames = ["sample", "segments", "stdev", "mad", "iqr", "bivar"]
    return pd.DataFrame.from_records(rows, columns=colnames)


def zip_repeater(iterable, repeatable):
    """Repeat a single segmentation to match the number of copy ratio inputs"""
    rpt_len = len(repeatable)
    if rpt_len == 1:
        rpt = repeatable[0]
        for it in iterable:
            yield it, rpt
    else:
        i = -1
        for i, (it, rpt) in enumerate(zip(iterable, repeatable)):
            yield it, rpt
        # Require lengths to match
        if i + 1 != rpt_len:
            raise ValueError("""Number of unsegmented and segmented input files
                             did not match (%d vs. %d)""" % (i, rpt_len))


def ests_of_scale(deviations):
    """Estimators of scale: standard deviation, MAD, biweight midvariance.

    Calculates all of these values for an array of deviations and returns them
    as a tuple.
    """
    std = np.std(deviations, dtype=np.float64)
    mad = descriptives.median_absolute_deviation(deviations)
    iqr = descriptives.interquartile_range(deviations)
    biw = descriptives.biweight_midvariance(deviations)
    return (std, mad, iqr, biw)
