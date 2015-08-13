"""NGS utilities: BED/Interval I/O, genomic regions and ranges."""
from __future__ import absolute_import, division, print_function

import functools
import re
import shlex

from Bio.File import as_handle

from .shared import echo


def sniff_num_columns(bed_fname):
    """Detect the number of columns in a BED/interval file.

    Guidance:
        3 cols => coordinates only;
        5 cols => intervals file (coordinates, strand, name);
        otherwise => Full or extended BED format
    """
    for firstrow in parse_regions(bed_fname):
        return len(firstrow)


def sniff_region_format(fname):
    """Guess whether the file format is BED, Picard interval list, or text.

    Returns a tuple of the format name (str) or None if the file is empty.
    """
    with open(fname, 'rU') as handle:
        for line in handle:
            if not line.strip():
                continue
            if '\t' not in line and ':' in line and '-' in line:
                return 'text'
            if line.startswith('@') or re.match('\w+\t\d+\t\d+\t(\+|-|\.)\t\S+',
                                                line):
                return 'interval'
            if line.startswith('track') or line.count('\t') > 1:
                return 'bed'
            raise ValueError("File " + repr(fname) + " does not appear to "
                             + "be BED, interval list, or 'chr:start-end' "
                             + "text!\nFirst non-blank line: " + repr(line))


def parse_regions(fname, coord_only=False, keep_strand=False):
    """Parse regions in any of the expected file formats.

    Iterates over tuples of the tabular contents. Header lines are skipped.

    Start and end coordinates are base-0, half-open.

    If coord_only, yield triplets of (chrom, start, end). Otherwise, yield
    quads of (chrom, start, end, name).
    """
    fmt = sniff_region_format(fname)
    if fmt is None:
        return []
    if fmt == 'bed':
        echo("Detected file format: BED")
    elif fmt == 'interval':
        echo("Detected file format: interval list")
    parser = {'text': parse_text_coords,
              'interval': parse_interval_list,
              'bed': parse_bed,
             }[fmt]
    return parser(fname, coord_only, keep_strand)


def report_bad_line(line_parser):
    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError:
            raise ValueError("Bad line: %r" % line)
    return wrapper


def parse_text_coords(fname, coord_only, _keep_strand):
    """Parse text coordinates: chrom:start-end

    Text coordinates are assumed to be counting from 1.
    """
    if coord_only:
        @report_bad_line
        def _parse_line(line):
            chrom, _rest = line.rstrip().split(':', 1)
            start, end = _rest.split('-')
            if ':' in end:
                end = end.split(':', 1)[0]
            return chrom, int(start) - 1, int(end)
    else:
        @report_bad_line
        def _parse_line(line):
            fields = line.split(':')
            if len(fields) == 3:
                chrom, start_end, name = fields
            elif len(fields) == 2:
                chrom, start_end = fields
                name = '-'
            else:
                raise ValueError
            start, end = start_end.split('-')
            return chrom, int(start) - 1, int(end), name.rstrip()

    with as_handle(fname, 'rU') as handle:
        for line in handle:
            yield _parse_line(line)


def parse_interval_list(fname, coord_only, keep_strand):
    """Parse a Picard-compatible interval list.

    Expected tabular columns:
        chromosome, start position, end position, strand, region name

    Counting is from 1.
    """
    if coord_only:
        if keep_strand:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end, strand = line.split('\t')[:4]
                return chrom, int(start) - 1, int(end), strand.rstrip()
        else:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end = line.split('\t')[:3]
                return chrom, int(start) - 1, int(end)
    elif keep_strand:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t')
            chrom, start, end, strand = fields[:4]
            if len(fields) > 4:
                name = fields[-1].rstrip()
            else:
                name = '-'
            return chrom, int(start) - 1, int(end), name, strand
    else:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t')
            chrom, start, end = fields[:3]
            if len(fields) > 3:
                name = fields[-1].rstrip()
            else:
                name = '-'
            return chrom, int(start) - 1, int(end), name

    with as_handle(fname, 'rU') as handle:
        for line in handle:
            if line.startswith('@'):
                # Skip the SAM header
                continue
            yield _parse_line(line)


def parse_bed(fname, coord_only, keep_strand):
    """Parse a BED file.

    A BED file has these columns:
        chromosome, start position, end position, [name, strand, other stuff...]

    Counting is from 0.

    Sets of regions are separated by "track" lines. This function stops
    iteration after encountering a track line other than the first one in the
    file.
    """
    if coord_only:
        if keep_strand:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end, _name, _score, strand = line.split('\t', 6)[:6]
                return chrom, int(start), int(end), strand.rstrip()
        else:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end = line.split('\t', 3)[:3]
                return chrom, int(start), int(end)
    elif keep_strand:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t', 6)
            chrom, start, end = fields[:3]
            name = (fields[3].rstrip()
                    if len(fields) >= 4 else '-')
            strand = (fields[5].rstrip()
                      if len(fields) >= 6 else '.')
            return chrom, int(start), int(end), name, strand
    else:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t', 4)
            chrom, start, end = fields[:3]
            name = (fields[3].rstrip()
                    if len(fields) >= 4 else '-')
            return chrom, int(start), int(end), name

    with as_handle(fname, 'rU') as handle:
        firstline = next(handle)
        if firstline.startswith("track"):
            pass
        else:
            yield _parse_line(firstline)

        for line in handle:
            if line.startswith('track'):
                raise StopIteration
            yield _parse_line(line)


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

