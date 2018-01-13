#!/usr/bin/env python

"""List the locations of accessible sequence regions in a FASTA file.

Inaccessible regions, e.g. telomeres and centromeres, are masked out with N in
the reference genome sequence; this script scans those to identify the
coordinates of the accessible regions (those between the long spans of N's).
"""
from __future__ import absolute_import, division, print_function
from builtins import next, zip

import logging

import numpy as np
from skgenome import tabio, GenomicArray as GA


def do_access(fa_fname, exclude_fnames=(), min_gap_size=5000,
              skip_noncanonical=True):
    """List the locations of accessible sequence regions in a FASTA file."""
    fa_regions = get_regions(fa_fname)
    if skip_noncanonical:
        fa_regions = drop_noncanonical_contigs(fa_regions)
    access_regions = GA.from_rows(fa_regions)
    for ex_fname in exclude_fnames:
        excluded = tabio.read(ex_fname, 'bed3')
        access_regions = access_regions.subtract(excluded)
    return GA.from_rows(join_regions(access_regions, min_gap_size))


def drop_noncanonical_contigs(region_tups):
    """Drop contigs with noncanonical names.

    `region_tups` is an iterable of (chrom, start, end) tuples.

    Yield the same, but dropping noncanonical chrom.
    """
    from .antitarget import is_canonical_contig_name
    return (tup for tup in region_tups
            if is_canonical_contig_name(tup[0]))


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
                logging.info("%s: Scanning for accessible regions", chrom)
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
                        n_indices = np.where(line_chars == b'N')[0]
                        # Emit the first block of non-N chars, if any
                        if run_start is not None:
                            yield log_this(chrom, run_start, cursor + n_indices[0])
                        elif n_indices[0] != 0:
                            yield log_this(chrom, cursor, cursor + n_indices[0])
                        # Emit any short intermediate blocks
                        gap_mask = np.diff(n_indices) > 1
                        if gap_mask.any():
                            ok_starts = n_indices[:-1][gap_mask] + 1 + cursor
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


def log_this(chrom, run_start, run_end):
    """Log a coordinate range, then return it as a tuple."""
    logging.info("\tAccessible region %s:%d-%d (size %d)",
                 chrom, run_start, run_end, run_end - run_start)
    return (chrom, run_start, run_end)


def join_regions(regions, min_gap_size):
    """Filter regions, joining those separated by small gaps."""
    min_gap_size = min_gap_size or 0
    for chrom, rows in regions.by_chromosome():
        logging.info("%s: Joining over small gaps", chrom)
        coords = iter(zip(rows['start'], rows['end']))
        prev_start, prev_end = next(coords)
        for start, end in coords:
            gap = start - prev_end
            assert gap > 0, ("Impossible gap between %s %d-%d and %d-%d (=%d)"
                             % (chrom, prev_start, prev_end, start, end, gap))
            if gap < min_gap_size:
                # Join with the previous region
                logging.info("\tJoining %s %d-%d and %d-%d (gap size %d)",
                             chrom, prev_start, prev_end, start, end, gap)
                prev_end = end
            else:
                # Keep the gap; emit the previous region as-is
                logging.info("\tKeeping gap %s:%d-%d (size %d)",
                             chrom, prev_end, start, gap)
                yield (chrom, prev_start, prev_end)
                prev_start, prev_end = start, end
        yield (chrom, prev_start, prev_end)
