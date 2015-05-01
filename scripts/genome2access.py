#!/usr/bin/env python

"""List the locations of accessible sequence regions in a FASTA file.

Inaccessible regions, e.g. telomeres and centromeres, are masked out with N in
the reference genome sequence; this script scans those to identify the
coordinates of the accessible regions (those between the long spans of N's).
"""

import argparse
import sys
import itertools

import numpy as np

from cnvlib.ngfrills import echo, parse_regions


def log_this(chrom, run_start, run_end):
    echo("\tAccessible region %s:%d-%d (size %d)"
         % (chrom, run_start, run_end, run_end - run_start))
    return (chrom, run_start, run_end)


def get_regions(fasta_fname):
    """Find accessible sequence regions (those not masked out with 'N')."""
    with open(fasta_fname) as infile:
        chrom = cursor = run_start = None
        for line in infile:
            if line.startswith('>'):
                # Emit the last chromosome's last run, if any
                if run_start is not None:
                    yield log_this(chrom, run_start, cursor)
                # Start new chromosome
                chrom = line.split(None, 1)[0][1:]
                run_start = None
                cursor = 0
                echo(chrom + ": Scanning for accessible regions")
            else:
                line = line.rstrip()
                if 'N' in line:
                    if all(c == 'N' for c in line):
                        # Shortcut if the line is all N chars
                        if run_start is not None:
                            yield log_this(chrom, run_start, cursor)
                            run_start = None
                    else:
                        # Slow route: line is a mix of N and non-N chars
                        line_chars = np.array(line, dtype='c')
                        n_indices = np.where(line_chars == 'N')[0]
                        # Emit the first block of non-N chars, if any
                        if run_start is not None:
                            yield log_this(chrom, run_start, cursor + n_indices[0])
                        elif n_indices[0] != 0:
                            yield log_this(chrom, cursor, cursor + n_indices[0])
                        # Emit any short intermediate blocks
                        gap_mask = np.diff(n_indices) > 1
                        if gap_mask.any():
                            ok_starts = n_indices[gap_mask] + 1 + cursor
                            ok_ends = n_indices[1:][gap_mask] + cursor
                            for start, end in zip(ok_starts, ok_ends):
                                yield log_this(chrom, start, end)
                        # Account for any tailing non-N chars
                        if n_indices[-1] + 1 < len(line_chars):
                            run_start = cursor + n_indices[-1] + 1
                        else:
                            run_start = None
                else:
                    if run_start is None:
                        # Start of a new run of non-N characters
                        run_start = cursor
                cursor += len(line)
        # Emit the last run if it's accessible (i.e. not a telomere)
        if run_start is not None:
            yield log_this(chrom, run_start, cursor)


def group_regions_by_chromosome(rows):
    """Iterate through BED3 rows: (chrom, BED3-rows-in-this-chrom)"""
    for chrom, rows in itertools.groupby(rows, lambda bed3: bed3[0]):
        yield chrom, [(start, end) for _chrom, start, end in rows]


def join_regions(regions, min_gap_size):
    """Filter regions, joining those separated by small gaps."""
    for chrom, coords in group_regions_by_chromosome(regions):
        echo(chrom + ": Joining over small gaps")
        coords = iter(coords)
        prev_start, prev_end = next(coords)
        for start, end in coords:
            gap = start - prev_end
            assert gap > 0, ("Impossible gap between %s %d-%d and %d-%d (=%d)"
                             % (chrom, prev_start, prev_end, start, end, gap))
            if gap < min_gap_size:
                # Join with the previous region
                echo("\tJoining %s %d-%d and %d-%d (gap size %d)"
                     % (chrom, prev_start, prev_end, start, end, gap))
                prev_end = end
            else:
                # Keep the gap; emit the previous region as-is
                echo("\tKeeping gap %s:%d-%d (size %d)"
                     % (chrom, prev_end, start, gap))
                yield (chrom, prev_start, prev_end)
                prev_start, prev_end = start, end
        yield (chrom, prev_start, prev_end)


def exclude_regions(bed_fname, access_rows):
    ex_by_chrom = dict(group_regions_by_chromosome(
        parse_regions(bed_fname, coord_only=True)))
    if len(ex_by_chrom) == 0:
        # Nothing to exclude -> emit the input regions unmodified
        for row in access_rows:
            yield row
    else:
        # Check if each input region overlaps an excluded region
        for chrom, a_rows in group_regions_by_chromosome(access_rows):
            if chrom in ex_by_chrom:
                echo(chrom + ": Subtracting excluded regions")
                exclude_rows = iter(ex_by_chrom[chrom])
                ex_start, ex_end = next_or_inf(exclude_rows)
                for a_start, a_end in a_rows:
                    for row in exclude_in_region(exclude_rows, chrom, a_start,
                                                 a_end, ex_start, ex_end):
                        yield row
            else:
                echo(chrom + ": No excluded regions")
                for a_start, a_end in a_rows:
                    yield (chrom, a_start, a_end)


def exclude_in_region(exclude_rows, chrom, a_start, a_end, ex_start, ex_end):
    """Take region exclusions from an iterable and apply, perhaps recursively.

    Returns an iterable (usually length 1) of two tuples:
        (accessible chromosome, start, end)
        (current exclusion start, end)
    """
    # If we've leapfrogged the excluded area, catch up
    while ex_end <= a_start:
        ex_start, ex_end = next_or_inf(exclude_rows)
    if a_end <= ex_start:
        # Excluded area does not overlap this one
        yield (chrom, a_start, a_end)
    else:
        # Excluded area overlaps this one -> trim this region
        echo("\tExclusion %s:%d-%d overlaps accessible region %d-%d"
             % (chrom, ex_start, ex_end, a_start, a_end))
        if ex_start <= a_start:
            if ex_end < a_end:
                # Exclusion covers this region's left (start) edge only
                for row in exclude_in_region(exclude_rows, chrom, ex_end, a_end,
                                             ex_start, ex_end):
                    yield row
            # Otherwise: Exclusion covers the whole region
        else:
            yield (chrom, a_start, ex_start)
            if ex_end < a_end:
                # Exclusion is in the middle of this region
                for row in exclude_in_region(exclude_rows, chrom, ex_end,
                                             a_end, ex_start, ex_end):
                    yield row
            # Otherwise: Exclusion covers this region's right (end) edge only


def next_or_inf(iterable):
    try:
        return next(iterable)
    except StopIteration:
        return float("Inf"), float("Inf")


if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("fa_fname",
                    help="Genome FASTA file name")
    AP.add_argument("-s", "--min-gap-size", type=int, default=5000,
                    help="""Minimum gap size between accessible sequence
                    regions.  Regions separated by less than this distance will
                    be joined together. [Default: %(default)s]""")
    AP.add_argument("-x", "--exclude", action="append", default=[],
                    help="""Additional regions to exclude, in BED format. Can be
                    used multiple times.""")
    AP.add_argument("-o", "--output",
                    type=argparse.FileType('w'), default=sys.stdout,
                    help="Output file name")
    args = AP.parse_args()

    # Closes over args.output
    def write_row(chrom, run_start, run_end):
        args.output.write("%s\t%s\t%s\n" % (chrom, run_start, run_end))
        args.output.flush()

    access_regions = get_regions(args.fa_fname)
    for ex_fname in args.exclude:
        access_regions = exclude_regions(ex_fname, access_regions)
    for row in join_regions(access_regions, args.min_gap_size):
        write_row(*row)

