"""Trivial segmentation: one segment per chromosome arm.

This procedure calculates a single segment's mean and coordinates for each
chromosome arm in the given array (splitting via CNA.by_arm()).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..segmetrics import segment_mean

if TYPE_CHECKING:
    from ..cnary import CopyNumArray


def segment_none(cnarr: CopyNumArray) -> CopyNumArray:
    """Return one trivial segment per chromosome arm in the given array."""
    colnames = ["chromosome", "start", "end", "log2", "gene", "probes"]
    rows = [
        (
            arm.chromosome.iat[0],
            arm.start.iat[0],
            arm.end.iat[-1],
            segment_mean(arm),
            "-",
            len(arm),
        )
        for _label, arm in cnarr.by_arm()
        if len(arm)
    ]
    table = pd.DataFrame.from_records(rows, columns=colnames)
    segarr = cnarr.as_dataframe(table)
    segarr.sort_columns()
    return segarr
