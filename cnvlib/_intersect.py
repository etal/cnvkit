"""DataFrame-level intersection operations.

Calculate overlapping regions, similar to bedtools intersect.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
from __future__ import print_function, absolute_import, division

import numpy as np


def _iter_ranges(table, chrom, starts, ends, mode):
    """Iterate through sub-ranges."""
    assert mode in ('inner', 'outer', 'trim')
    # Optional if we've already subsetted by chromosome (not checked!)
    if chrom:
        assert isinstance(chrom, basestring)  # ENH: accept array?
        try:
            table = table[table['chromosome'] == chrom]
        except KeyError:
            raise KeyError("Chromosome %s is not in this probe set" % chrom)

    # Edge cases
    if not len(table):
        yield []
        raise StopIteration

    if starts is None and ends is None:
        yield table
        raise StopIteration

    if starts is not None and len(starts):
        if mode == 'inner':
            # Only rows entirely after the start point
            start_idxs = table.start.searchsorted(starts)
        else:
            # Include all rows overlapping the start point
            start_idxs = table.end.searchsorted(starts, 'right')
    else:
        starts = np.zeros(len(ends) if ends is not None else 1,
                            dtype=np.int_)
        start_idxs = starts.copy()

    if ends is not None and len(ends):
        if mode == 'inner':
            end_idxs = table.end.searchsorted(ends, 'right')
        else:
            end_idxs = table.start.searchsorted(ends)
    else:
        end_idxs = np.repeat(len(table), len(starts))
        ends = [None] * len(starts)

    for start_idx, start_val, end_idx, end_val in zip(start_idxs, starts,
                                                        end_idxs, ends):
        subtable = table[start_idx:end_idx]
        if mode == 'trim':
            subtable = subtable.copy()
            # Update 5' endpoints to the boundary
            if start_val:
                subtable.start = subtable.start.clip_lower(start_val)
            # Update 3' endpoints to the boundary
            if end_val:
                subtable.end = subtable.end.clip_upper(end_val)
        yield subtable
