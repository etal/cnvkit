"""I/O for UCSC Browser Extensible Data (BED)."""
from __future__ import absolute_import, division, print_function
from builtins import map, next

import shlex

import pandas as pd
from Bio.File import as_handle

from .util import report_bad_line


def read_bed(infile):
    """UCSC Browser Extensible Data (BED) format.

    A BED file has these columns:
        chromosome, start position, end position, [gene, strand, other stuff...]

    Coordinate indexing is from 0.

    Sets of regions are separated by "track" lines. This function stops reading
    after encountering a track line other than the first one in the file.
    """
    # ENH: just pd.read_table, skip 'track'
    @report_bad_line
    def _parse_line(line):
        fields = line.split('\t', 6)
        chrom, start, end = fields[:3]
        gene = (fields[3].rstrip()
                if len(fields) >= 4 else '-')
        strand = (fields[5].rstrip()
                if len(fields) >= 6 else '.')
        return chrom, int(start), int(end), gene, strand

    def track2track(handle):
        try:
            firstline = next(handle)
        except StopIteration:
            pass
        else:
            if not firstline.startswith("track"):
                yield firstline
            for line in handle:
                if line.startswith("track"):
                    break
                yield line

    with as_handle(infile, 'rU') as handle:
        rows = map(_parse_line, track2track(handle))
        return pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                        "end", "gene", "strand"])


def read_bed3(infile):
    """3-column BED format: chromosome, start, end."""
    table = read_bed(infile)
    return table.loc[:, ['chromosome', 'start', 'end']]


def read_bed4(infile):
    """4-column BED format: chromosome, start, end, name."""
    table = read_bed(infile)
    return table.loc[:, ['chromosome', 'start', 'end', 'gene']]


def read_bed6(infile):
    """6-column BED format: chromosome, start, end, name, score, strand."""
    return NotImplemented


def parse_bed_track(line):
    """Parse the "name" field of a BED track definition line.

    Example:
    track name=146793_BastianLabv2_P2_target_region description="146793_BastianLabv2_P2_target_region"
    """
    fields = shlex.split(line)  # raises ValueError if line is corrupted
    assert fields[0] == 'track'
    for field in fields[1:]:
        if '=' in field:
            key, val = field.split('=', 1)
            if key == 'name':
                return val
    raise ValueError("No name defined for this track")


def group_bed_tracks(bedfile):
    """Group the parsed rows in a BED file by track.

    Yields (track_name, iterable_of_lines), much like itertools.groupby.
    """
    # ENH - make this memory-efficient w/ generators or something
    with as_handle(bedfile, 'r') as handle:
        curr_track = 'DEFAULT'
        curr_lines = []
        for line in handle:
            if line.startswith('track'):
                if curr_lines:
                    yield curr_track, curr_lines
                    curr_lines = []
                curr_track = parse_bed_track(line)
            else:
                curr_lines.append(line)
        yield curr_track, curr_lines


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
