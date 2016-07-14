"""I/O for tab-delimited format with an initial header row.

Column names match the target class attributes. At least "chromosome", "start",
and "end" are required.
"""
from __future__ import absolute_import, division, print_function

import pandas as pd


def read_tab(infile):
    """Read tab-separated data with column names in the first row.

    The format is BED-like, but with a header row included and with
    arbitrary extra columns.
    """
    return pd.read_table(infile, dtype={'chromosome': 'str'})


def write_tab(dframe):
    """Write tab-separated data with column names in the first row."""
    return dframe

