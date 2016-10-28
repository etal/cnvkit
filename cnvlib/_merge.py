"""DataFrame-level merging operations.

Merge overlapping regions into single rows, similar to bedtools merge.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
from __future__ import print_function, absolute_import, division

import numpy as np
import pandas as pd

    # NB: same gene name can appear on alt. contigs

# Merge overlapping rows
# XXX see also target.py

def _merge(table, stranded=False, combiners=None):
    cmb = {
        'start': lambda ser: ser.iat[0], # pd.Series.min,
        'end': lambda ser: ser.iat[-1], # pd.Series.max,
        'gene': join_strings,
        'accession': join_strings,
    }
    if combiners:
        cmb.update(combiners)
    if stranded:
        groupkey = ['chromosome', 'strand']
        if 'strand' not in combiners:
            combiners['strand'] = lambda ser: ser.iat[0]
    else:
        groupkey = ['chromosome']
        if 'strand' not in combiners:
            combiners['strand'] = merge_strands
    return (table.groupby(by=groupkey,
                          as_index=False, group_keys=False, sort=False)
            .apply(_merge_overlapping, cmb))


def _merge_overlapping(table, combiners):
    """Merge overlapping regions within a group."""
    # Short-circuit the simple, common cases
    if len(table) == 1:
        return table
    if table['start'].nunique() == table['end'].nunique() == 1:
        return _squash_rows(table, combiners)
    # Identify & enumerate (non)overlapping groups of rows
    overlap_sizes = table.end.cummax().values[:-1] - table.start.values[1:]
    non_overlapping = np.r_[False, (overlap_sizes <= 0)]
    # Squash rows within each non-overlapping group
    return (table.groupby(non_overlapping * np.arange(len(non_overlapping)),
                          as_index=False, group_keys=False, sort=False)
            .apply(_squash_rows, combiners))


# Squash rows according to a given grouping criterion
# XXX see also segfilter.py

def _squash_rows(table, combiners):
    """Reduce multiple rows into one, combining 'accession' field."""
    i = table.first_valid_index()
    row = table.loc[i:i, :]
    for key, combiner in combiners.viewitems():
        row[key] = combiner(table[key])
    return row


# Combiners

def join_strings(ser):
    """Join a Series of strings by commas."""
    # TODO if ser elements are also comma-separated, split+uniq those too
    return ','.join(ser.unique())


def merge_strands(ser):
    strands = ser.drop_duplicates()
    if len(strands) > 1:
        return '.'
    return strands.iat[0]
