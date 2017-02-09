from __future__ import absolute_import, division, print_function
import collections

import pandas as pd
from Bio.File import as_handle

from .util import report_bad_line
from ..genome import GenomicArray as GA


def read_text(infile):
    """Text coordinate format: "chr:start-end", one per line.

    Or sometimes: "chrom:start-end gene" or "chrom:start-end REF>ALT"

    Coordinate indexing is assumed to be from 1.
    """
    parse_line = report_bad_line(from_label)
    with as_handle(infile, 'rU') as handle:
        rows = [parse_line(line) for line in handle]
    return pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                    "end", "gene"])


Region = collections.namedtuple('Region', 'chromosome start end gene')

def from_label(label):
    fields = label.split(':')
    if len(fields) == 3:
        chrom, start_end, gene = fields
    elif len(fields) == 2:
        chrom, start_end = fields
        gene = '-'
    else:
        raise ValueError("Bad label: %r" % label)
    start, end = start_end.split('-')
    return Region(chrom, int(start) - 1, int(end), gene.rstrip())

# ENH: import from_label, to_label at top level in tabio.__init__?
# - move GA.row2label implementation here? Both import it from core/util?
to_label = GA.row2label


def write_text(dframe):
    """Text coordinate format: "chr:start-end", one per line."""
    dframe = dframe.copy()
    dframe['start'] += 1
    return dframe.apply(GA.row2label, axis=1)
