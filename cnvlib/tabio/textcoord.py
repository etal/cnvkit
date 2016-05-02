from __future__ import absolute_import, division, print_function

import pandas as pd
from Bio.File import as_handle

from .. import ngfrills
from ..gary import GenomicArray as GA


def read_text(infile):
    """Text coordinate format: "chr:start-end", one per line.

    Or sometimes: "chrom:start-end gene" or "chrom:start-end REF>ALT"

    Coordinate indexing is assumed to be from 1.
    """
    @ngfrills.report_bad_line
    def _parse_line(line):
        fields = line.split(':')
        if len(fields) == 3:
            chrom, start_end, gene = fields
        elif len(fields) == 2:
            chrom, start_end = fields
            gene = '-'
        else:
            raise ValueError("Bad line: %r" % line)
        start, end = start_end.split('-')
        return chrom, int(start) - 1, int(end), gene.rstrip()

    with as_handle(infile, 'rU') as handle:
        rows = [_parse_line(line) for line in handle]
    return pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                    "end", "gene"])


def write_text(dframe):
    """Text coordinate format: "chr:start-end", one per line."""
    dframe = dframe.copy()
    dframe['start'] += 1
    return dframe.apply(GA.row2label, axis=1)
