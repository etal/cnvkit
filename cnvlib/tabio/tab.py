"""I/O for tab-delimited format with an initial header row.

Column names match the target class attributes. At least "chromosome", "start",
and "end" are required.
"""
from __future__ import absolute_import, division, print_function

import logging

import pandas as pd


def read_tab(infile):
    """Read tab-separated data with column names in the first row.

    The format is BED-like, but with a header row included and with
    arbitrary extra columns.
    """
    dframe = pd.read_table(infile, dtype={'chromosome': 'str'})
    if "log2" in dframe.columns:
        # Every bin needs a log2 value; the others can be NaN
        d2 = dframe.dropna(subset=["log2"])
        if len(d2) < len(dframe):
            logging.warn("Dropped %d rows with missing log2 values",
                        len(dframe) - len(d2))
            dframe = d2.copy()
    return dframe


def write_tab(dframe):
    """Write tab-separated data with column names in the first row."""
    return dframe

