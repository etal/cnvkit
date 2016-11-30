"""DataFrame-level intersection operations.

Calculate overlapping regions, similar to bedtools intersect.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
from __future__ import print_function, absolute_import, division
from builtins import zip
from past.builtins import basestring

import logging

import numpy as np
import pandas as pd

from .combiners import first_of, join_strings, make_const


def by_ranges(table, other, mode, keep_empty):
    """Group rows by another GenomicArray's bin coordinate ranges."""
    for _chrom, bin_rows, src_rows in by_shared_chroms(other, table,
                                                       keep_empty):
        if src_rows is not None:
            subranges = iter_ranges(src_rows, None, bin_rows['start'],
                                    bin_rows['end'], mode)
            for bin_row, subrange in zip(bin_rows.itertuples(index=False),
                                         subranges):
                yield bin_row, subrange
        elif keep_empty:
            for bin_row in bin_rows.itertuples(index=False):
                yield bin_row, []  # ENH: empty dframe matching table


def by_shared_chroms(table, other, keep_empty=True):
    other_chroms = {c: o for c, o in other.groupby(['chromosome'], sort=False)}
    for chrom, ctable in table.groupby(['chromosome'], sort=False):
        if chrom in other_chroms:
            otable = other_chroms[chrom]
            yield chrom, ctable, otable
        elif keep_empty:
            yield chrom, ctable, None


def into_ranges(source, dest, src_col, default, summary_func):
    """Assign to `dest` from aligned rows in `source`."""
    if summary_func is None:
        # Choose a type-appropriate summary function
        elem = source[src_col].iat[0]
        if isinstance(elem, (basestring, np.string_)):
            summary_func = join_strings
        elif isinstance(elem, (float, np.float_)):
            summary_func = np.nanmedian
        else:
            summary_func = first_of
    elif not callable(summary_func):
        # Just fill in the given value, I suppose.
        summary_func = make_const(summary_func)

    def series2value(ser):
        if len(ser) == 1:
            return ser.iat[0]
        return summary_func(ser)

    return pd.Series([(series2value(src_rows[src_col])
                       if len(src_rows) else default)
                      for _bin, src_rows in by_ranges(source, dest, 'outer',
                                                      True)])


def iter_ranges(table, chrom, starts, ends, mode):
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
