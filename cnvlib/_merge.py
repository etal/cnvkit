"""DataFrame-level merging operations.

Merge overlapping regions into single rows, similar to bedtools merge.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
from __future__ import print_function, absolute_import, division

import itertools

import numpy as np
import pandas as pd


# Merge overlapping rows
# XXX see also target.py

def _merge(table, stranded=False, combiners=None):
    cmb = {
        # # For pandas.Series instances
        # 'start': lambda ser: ser.iat[0],
        # 'end': lambda ser: ser.iat[-1],
        # For Python lists
        'start': lambda a: a[0],
        'end': lambda a: a[-1],
        'gene': join_strings,
        'accession': join_strings,
    }
    if combiners:
        cmb.update(combiners)
    if stranded:
        groupkey = ['chromosome', 'strand']
        if 'strand' not in cmb:
            # cmb['strand'] = lambda ser: ser.iat[0]  # pd.Series
            cmb['strand'] = lambda a: a[0]  # Py list
    else:
        # NB: same gene name can appear on alt. contigs
        groupkey = ['chromosome']
        if 'strand' not in cmb:
            cmb['strand'] = merge_strands
    table = table.sort_values(groupkey + ['start', 'end'])
    return (table.groupby(by=groupkey,
                          as_index=False, group_keys=False, sort=False)
            .apply(_merge_overlapping, cmb))


def _merge_overlapping(table, combiners):
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
    keyed_groups = itertools.izip(_nonoverlapping_groups(table),
                                  table.itertuples(index=False))
    merged_rows = [_squash_tuples(row_group, combiners)
                   for _key, row_group in itertools.groupby(keyed_groups,
                                                            lambda x: x[0])]
    return pd.DataFrame.from_records(merged_rows, columns=merged_rows[0]._fields)


def _nonoverlapping_groups(table):
    """Identify and enumerate groups of overlapping rows.

    That is, increment the group ID after each non-negative gap between
    intervals. Intervals (rows) will be merged if any bases overlap.
    """
    # Examples:
    #
    #  gap?     F  T  F  F  T  T  T  F
    #  group  0  1  1  2  3  3  3  3  4
    #
    #  gap?     T  F  T  T  F  F  F  T
    #  group  0  0  1  1  1  2  3  4  4
    gap_sizes = table.start.values[1:] - table.end.cummax().values[:-1]
    return np.r_[False, gap_sizes >= 0].cumsum()


# Squash rows according to a given grouping criterion
# XXX see also segfilter.py

def _squash_tuples(keyed_rows, combiners):
    """Combine multiple rows into one NamedTuple."""
    rows = [kr[1] for kr in keyed_rows] #list(rows)
    firsttup = rows[0]
    if len(rows) == 1:
        return firsttup
    newfields = {key: combiner([getattr(r, key) for r in rows])
                 for key, combiner in combiners.viewitems()}
    return firsttup._replace(**newfields)


# Combiners

def join_strings(elems):
    """Join a Series of strings by commas."""
    # TODO if ser elements are also comma-separated, split+uniq those too
    return ','.join(pd.unique(elems))


def merge_strands(elems):
    strands = set(elems)
    if len(strands) > 1:
        return '.'
    return elems[0]


# Combiners -- for pandas.Series instances

# def join_strings(ser):
#     """Join a Series of strings by commas."""
#     # TODO if ser elements are also comma-separated, split+uniq those too
#     return ','.join(ser.unique())


# def merge_strands(ser):
#     strands = ser.drop_duplicates()
#     if len(strands) > 1:
#         return '.'
#     return strands.iat[0]
