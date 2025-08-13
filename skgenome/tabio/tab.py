"""I/O for tab-delimited format with an initial header row.

Column names match the target class attributes. At least "chromosome", "start",
and "end" are required.
"""

from __future__ import annotations
import logging

import pandas as pd
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pandas.core.frame import DataFrame


def read_tab(infile: str) -> DataFrame:
    """Read tab-separated data with column names in the first row.

    The format is BED-like, but with a header row included and with
    arbitrary extra columns.
    """
    dframe = pd.read_csv(infile, sep="\t", dtype={"chromosome": "str"})
    if "log2" in dframe.columns:
        # Every bin needs a log2 value; the others can be NaN
        d2 = dframe.dropna(subset=["log2"])
        if len(d2) < len(dframe):
            logging.warning(
                "Dropped %d rows with missing log2 values", len(dframe) - len(d2)
            )
            dframe = d2.copy()
    return dframe


def write_tab(dframe: DataFrame) -> DataFrame:
    """Write tab-separated data with column names in the first row."""
    return dframe
