"""I/O for formats used by Picard tools.

- Interval list (also used in GATK)
- CalculateHsMetrics PER_TARGET_COVERAGE output

"""
from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd


def read_interval(infile):
    """GATK/Picard-compatible interval list format.

    Expected tabular columns:
        chromosome, start position, end position, strand, gene

    Coordinate indexing is from 1.
    """
    dframe = pd.read_table(infile,
                           comment='@', # Skip the SAM header
                           names=["chromosome", "start", "end", "strand", "gene",
                                 ])
    dframe["gene"].fillna('-', inplace=True)
    dframe["start"] -= 1
    return dframe


def read_picard_hs(infile):
    """Picard CalculateHsMetrics PER_TARGET_COVERAGE.

    The format is BED-like, but with a header row and the columns::

        chrom (str),
        start, end, length (int),
        name (str),
        %gc, mean_coverage, normalized_coverage (float)

    """
    dframe = pd.read_table(infile, na_filter=False, dtype={
        "chrom": "str",
        "start": "int",
        "end": "int",
        "length": "int",
        "name": "str",
        "%gc": "float",
        "mean_coverage": "float",
        "normalized_coverage": "float",
    })
    dframe.columns = ["chromosome", # chrom
                      "start", "end", "length",
                      "gene", # name
                      "gc", # %gc
                      "depth", "ratio"]
    del dframe["length"]
    dframe["start"] -= 1
    return dframe


# _____________________________________________________________________

def write_interval(dframe):
    dframe = dframe.copy()
    dframe["start"] += 1
    if "gene" not in dframe:
        dframe["gene"] = '-'
    if "strand" not in dframe:
        dframe["strand"] = "+"
    return dframe.loc[:, ["chromosome", "start", "end", "strand", "gene"]]


def write_picard_hs(dframe):
    if "depth" in dframe.columns:
        coverage = dframe["depth"]
        norm = coverage / coverage.mean()
    else:
        coverage = np.exp2(dframe["log2"])
        norm = coverage
    return pd.DataFrame.from_items([
        ("chrom", dframe["chromosome"]),
        ("start", dframe["start"] + 1),
        ("end", dframe["end"]),
        ("length", dframe["end"] - dframe["start"]),
        ("name", dframe["gene"]),
        ("%gc", dframe["gc"]),
        ("mean_coverage", coverage),
        ("normalized_coverage", norm),
    ])

