"""Simple VCF I/O.

Read only coordinate info & store the remaining columns as unparsed strings.
Just enough functionality to extract a subset of samples and/or perform
bedtools-like operations on VCF records.
"""
from __future__ import absolute_import, division, print_function

import logging

import pandas as pd

from ..vary import VariantArray as VA

# TODO
#   save VCF header (as string, the whole text block) in meta{header=}
#   then splode the INFO column
#       1st     just as strings?
#       Then    parse according to VCF 4.3 spec (floats, tuples, flags, etc.)
#               matching pysam as much as possible
def read_vcf(infile, skip_reject=False):
    """Read VCF file w/o samples."""
    # ENH: use pysam reader to parse INFO column?
    columns = ['chromosome', 'start', 'ref', 'alt', # 'filter', 'info',
              ]
    dtypes = [str, int, str, str, # str, str
             ]
    table = pd.read_table(infile,
                          comment="#",
                          header=None,
                          na_filter=False,
                          names=["chromosome", "start", "_ID", "ref", "alt",
                                 "_QUAL", "filter", "info"],
                          usecols=columns,
                          # ENH: converters={'info': func to parse it}
                          dtype=dict(zip(columns, dtypes)),
                         )
    # ENH: do things with filter, info
    # if skip_reject and record.FILTER and len(record.FILTER) > 0:
    table['end'] = table['start'] + table["alt"].str.len()  # ENH: INFO["END"]
    table['start'] -= 1
    logging.info("Loaded %d plain records", len(table))
    return table.loc[:, VA._required_columns]
