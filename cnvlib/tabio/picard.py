"""I/O for formats used by Picard tools.

- Interval list (also used in GATK)
- CalculateHsMetrics PER_TARGET_COVERAGE output

"""
from __future__ import absolute_import, division, print_function
# from builtins import map, next

import logging

import numpy as np
import pandas as pd

from .. import params


TOO_MANY_NO_COVERAGE = 100

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

    The format is BED-like, but with a header row and the columns:

    -
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
    dframe["gene"] = dframe["gene"].apply(unpipe_name)
    # Avoid math domain error converting coverages to log2 scale
    coverages = dframe["ratio"].copy()
    no_cvg_idx = (coverages == 0)
    if no_cvg_idx.sum() > TOO_MANY_NO_COVERAGE:
        logging.warn("*WARNING* Sample %s has >%d bins with no coverage",
                     str(infile), TOO_MANY_NO_COVERAGE)
    coverages[no_cvg_idx] = 2**params.NULL_LOG2_COVERAGE
    dframe["log2"] = np.log2(coverages)
    return dframe


def unpipe_name(name):
    """Fix the duplicated gene names Picard spits out.

    Return a string containing the single gene name, sans duplications and pipe
    characters.

    Picard CalculateHsMetrics combines the labels of overlapping intervals
    by joining all labels with '|', e.g. 'BRAF|BRAF' -- no two distinct
    targeted genes actually overlap, though, so these dupes are redundant.
    Meaningless target names are dropped, e.g. 'CGH|FOO|-' resolves as 'FOO'.
    In case of ambiguity, the longest name is taken, e.g. "TERT|TERT Promoter"
    resolves as "TERT Promoter".
    """
    if '|' not in name:
        return name
    gene_names = set(name.split('|'))
    if len(gene_names) == 1:
        return gene_names.pop()
    cleaned_names = gene_names.difference(params.IGNORE_GENE_NAMES)
    if cleaned_names:
        gene_names = cleaned_names
    new_name = sorted(gene_names, key=len, reverse=True)[0]
    if len(gene_names) > 1:
        logging.warn("*WARNING* Ambiguous gene name %r; using %r",
                     name, new_name)
    return new_name



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

