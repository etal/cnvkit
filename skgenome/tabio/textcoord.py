from __future__ import absolute_import, division, print_function

import pandas as pd
from Bio.File import as_handle

from skgenome.rangelabel import from_label, to_label
from .util import report_bad_line


def read_text(infile):
    """Text coordinate format: "chr:start-end", one per line.

    Or sometimes: "chrom:start-end gene" or "chrom:start-end REF>ALT"

    Coordinate indexing is assumed to be from 1.
    """
    parse_line = report_bad_line(from_label)
    with as_handle(infile, 'rU') as handle:
        rows = [parse_line(line) for line in handle]
    table = pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                     "end", "gene"])
    table['gene'] = table['gene'].replace('', '-')
    return table


def write_text(dframe):
    """Text coordinate format: "chr:start-end", one per line."""
    dframe = dframe.copy()
    dframe['start'] += 1
    return dframe.apply(to_label, axis=1)
