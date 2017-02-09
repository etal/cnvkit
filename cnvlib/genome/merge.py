"""DataFrame-level merging operations.

Merge overlapping regions into single rows, similar to bedtools merge.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
from __future__ import print_function, absolute_import, division
from builtins import zip

import itertools

import numpy as np
import pandas as pd

from .chromsort import sorter_chrom
from .combiners import get_combiners, first_of


def flatten(table, combine=None):
    if not len(table):
        return table
    if (table.start.values[1:] >= table.end.cummax().values[:-1]).all():
        return table
    # NB: Input rows and columns should already be sorted like this
    table = table.sort_values(['chromosome', 'start', 'end'])
    cmb = get_combiners(table, False, combine)
    out = (table.groupby(by='chromosome',
                         as_index=False, group_keys=False, sort=False)
           .apply(_flatten_overlapping, cmb)
           .reset_index(drop=True))
    return out.reindex(out.chromosome.apply(sorter_chrom)
                       .sort_values(kind='mergesort').index)


def _flatten_overlapping(table, combine):
    """Merge overlapping regions within a chromosome/strand.

    Assume chromosome and (if relevant) strand are already identical, so only
    start and end coordinates are considered.
    """
    keyed_groups = zip(_nonoverlapping_groups(table, 0),
                       table.itertuples(index=False))
    #  flat_rows = itertools.chain(
    #      *(_flatten_tuples(row_group, combine)
    #        for _key, row_group in itertools.groupby(keyed_groups, first_of)))
    #  out =  pd.DataFrame.from_records(list(flat_rows),
    #                                   columns=table.columns)
    #                                   #  columns=flat_rows[0]._fields)
    bits = [pd.DataFrame.from_records(list(_flatten_tuples(row_group, combine)),
                                      columns=table.columns)
            for _key, row_group in itertools.groupby(keyed_groups, first_of)]
    out = pd.concat(bits)

    return out


def _flatten_tuples(keyed_rows, combine):
    """Divide multiple rows where they overlap.

    Parameters
    ----------
    keyed_rows : iterable
        pairs of (non-overlapping-group index, overlapping rows)
    combine : dict

    Returns
    -------
    DataFrame
    """
    # TODO speed this up!
    rows = [kr[1] for kr in keyed_rows]
    first_row = rows[0]
    if len(rows) == 1:
        yield first_row
        raise StopIteration

    extra_cols = first_row._fields[3:]
    breaks = sorted(set(itertools.chain(*[(r.start, r.end) for r in rows])))
    for bp_start, bp_end in zip(breaks[:-1], breaks[1:]):
        # Find the row(s) overlapping this region
        # i.e. any not already seen and not already passed
        rows_in_play = []
        for row in rows:
            if row.start <= bp_start and row.end >= bp_end:
                rows_in_play.append(row)

        # Combine the extra fields of the overlapping regions
        extra_fields = {key: combine[key]([getattr(r, key)
                                           for r in rows_in_play])
                        for key in extra_cols}
        yield first_row._replace(start=bp_start, end=bp_end,
                                 **extra_fields)


def merge(table, bp=0, stranded=False, combine=None):
    """Merge overlapping rows in a DataFrame."""
    if not len(table):
        return table
    gap_sizes = table.start.values[1:] - table.end.cummax().values[:-1]
    if (gap_sizes > -bp).all():
        return table

    if stranded:
        groupkey = ['chromosome', 'strand']
    else:
        # NB: same gene name can appear on alt. contigs
        groupkey = ['chromosome']
    table = table.sort_values(groupkey + ['start', 'end'])
    cmb = get_combiners(table, stranded, combine)
    out = (table.groupby(by=groupkey,
                         as_index=False, group_keys=False, sort=False)
           .apply(_merge_overlapping, bp, cmb)
           .reset_index(drop=True))
    # Re-sort chromosomes cleverly instead of lexicographically
    return out.reindex(out.chromosome.apply(sorter_chrom)
                       .sort_values(kind='mergesort').index)


def _merge_overlapping(table, bp, combine):
    """Merge overlapping regions within a chromosome/strand.

    Assume chromosome and (if relevant) strand are already identical, so only
    start and end coordinates are considered.
    """
    # NB: pandas groupby seems like the obvious choice over itertools, but it is
    # very slow -- probably because it traverses the whole table (i.e.
    # chromosome) again to select groups, redoing the various inferences that
    # would be worthwhile outside this inner loop. With itertools, we take
    # advantage of the grouping and sorting already done, and don't repeat
    # pandas' traversal and inferences.
    # ENH: Find & use a lower-level, 1-pass pandas function
    keyed_groups = zip(_nonoverlapping_groups(table, bp),
                       table.itertuples(index=False))
    merged_rows = [_squash_tuples(row_group, combine)
                   for _key, row_group in itertools.groupby(keyed_groups,
                                                            first_of)]
    return pd.DataFrame.from_records(merged_rows,
                                     columns=merged_rows[0]._fields)


def _nonoverlapping_groups(table, bp):
    """Identify and enumerate groups of overlapping rows.

    That is, increment the group ID after each non-negative gap between
    intervals. Intervals (rows) will be merged if any bases overlap by at least
    `bp`.
    """
    # Examples:
    #
    #  gap?     F  T  F  F  T  T  T  F
    #  group  0  1  1  2  3  3  3  3  4
    #
    #  gap?     T  F  T  T  F  F  F  T
    #  group  0  0  1  1  1  2  3  4  4
    gap_sizes = table.start.values[1:] - table.end.cummax().values[:-1]
    return np.r_[False, gap_sizes > (-bp)].cumsum()


# Squash rows according to a given grouping criterion
# XXX see also segfilter.py
def _squash_tuples(keyed_rows, combine):
    """Combine multiple rows into one NamedTuple.

    Parameters
    ----------
    keyed_rows : iterable
        pairs of (non-overlapping-group index, overlapping rows)
    combine : dict

    Returns
    -------
    namedtuple
    """
    rows = [kr[1] for kr in keyed_rows] #list(rows)
    firsttup = rows[0]
    if len(rows) == 1:
        return firsttup
    newfields = {key: combiner([getattr(r, key) for r in rows])
                 for key, combiner in combine.items()}
    return firsttup._replace(**newfields)
