from __future__ import absolute_import, division, print_function

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

# TODO - picard_hs_pertarget (Picard CalculateHsMetrics per_target_coverages)

# _____________________________________________________________________

def write_interval(dframe):
    dframe = dframe.copy()
    dframe["start"] += 1
    if "gene" not in dframe:
        dframe["gene"] = '-'
    if "strand" not in dframe:
        dframe["strand"] = "+"
    return dframe.loc[:, ["chromosome", "start", "end", "strand", "gene"]]

