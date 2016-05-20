"""I/O for UCSC Browser Extensible Data (BED)."""
from __future__ import absolute_import, division, print_function
from builtins import map, next

import pandas as pd
from Bio.File import as_handle

from .. import ngfrills


def read_bed(infile):
    """UCSC Browser Extensible Data (BED) format.

    A BED file has these columns:
        chromosome, start position, end position, [gene, strand, other stuff...]

    Coordinate indexing is from 0.

    Sets of regions are separated by "track" lines. This function stops reading
    after encountering a track line other than the first one in the file.
    """
    # ENH: just pd.read_table, skip 'track'
    @ngfrills.report_bad_line
    def _parse_line(line):
        fields = line.split('\t', 6)
        chrom, start, end = fields[:3]
        gene = (fields[3].rstrip()
                if len(fields) >= 4 else '-')
        strand = (fields[5].rstrip()
                if len(fields) >= 6 else '.')
        return chrom, int(start), int(end), gene, strand

    def track2track(handle):
        firstline = next(handle)
        if firstline.startswith("track"):
            pass
        else:
            yield firstline
        for line in handle:
            if line.startswith('track'):
                raise StopIteration
            yield line

    with as_handle(infile, 'rU') as handle:
        rows = map(_parse_line, track2track(handle))
        return pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                        "end", "gene", "strand"])


# TODO - would these be useful?
def read_bed3(infile):
    """3-column BED format: chromosome, start, end."""
    return NotImplemented

def read_bed4(infile):
    """4-column BED format: chromosome, start, end, name."""
    return NotImplemented

def read_bed6(infile):
    """6-column BED format: chromosome, start, end, name, score, strand."""
    return NotImplemented


# _____________________________________________________________________

def write_bed(dframe):
    if len(dframe.columns) == 3:
        return write_bed3(dframe)
    elif len(dframe.columns) == 3:
        return write_bed4(dframe)
    else:
        # Default: BED-like, keep all trailing columns
        return dframe


def write_bed3(dframe):
    return dframe.loc[:, ["chromosome", "start", "end"]]


def write_bed4(dframe):
    dframe = dframe.copy()
    if "gene" not in dframe:
        dframe["gene"] = '-'
    return dframe.loc[:, ["chromosome", "start", "end", "gene"]]

