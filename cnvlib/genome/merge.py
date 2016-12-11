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
from .combiners import first_of, last_of, join_strings, merge_strands


def merge(table, stranded=False, combine=None):
    """Merge overlapping rows in a DataFrame."""
    cmb = {
        'start': first_of,
        'end': last_of,
        'gene': join_strings,
        'accession': join_strings,
    }
    if combine:
        cmb.update(combine)
    if stranded:
        groupkey = ['chromosome', 'strand']
        if 'strand' not in cmb:
            cmb['strand'] = first_of
    else:
        # NB: same gene name can appear on alt. contigs
        groupkey = ['chromosome']
        if 'strand' not in cmb:
            cmb['strand'] = merge_strands
    table = table.sort_values(groupkey + ['start', 'end'])
    cmb = {k: v for k, v in cmb.items() if k in table}
    out = (table.groupby(by=groupkey,
                         as_index=False, group_keys=False, sort=False)
           .apply(_merge_overlapping, cmb)
           .reset_index(drop=True))
    # Re-sort chromosomes cleverly instead of lexicographically
    return out.reindex(out.chromosome.apply(sorter_chrom)
                       .sort_values(kind='mergesort').index)


def _merge_overlapping(table, combine):
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
    keyed_groups = zip(_nonoverlapping_groups(table),
                       table.itertuples(index=False))
    merged_rows = [_squash_tuples(row_group, combine)
                   for _key, row_group in itertools.groupby(keyed_groups,
                                                            first_of)]
    return pd.DataFrame.from_records(merged_rows,
                                     columns=merged_rows[0]._fields)


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

def _squash_tuples(keyed_rows, combine):
    """Combine multiple rows into one NamedTuple."""
    rows = [kr[1] for kr in keyed_rows] #list(rows)
    firsttup = rows[0]
    if len(rows) == 1:
        return firsttup
    newfields = {key: combiner([getattr(r, key) for r in rows])
                 for key, combiner in combine.items()}
    return firsttup._replace(**newfields)
