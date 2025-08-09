"""DataFrame-level merging operations.

Merge overlapping regions into single rows, similar to bedtools merge.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
import itertools
from typing import Callable, Dict, Optional
from collections.abc import Iterable

import numpy as np
import pandas as pd

from .chromsort import sorter_chrom
from .combiners import get_combiners, first_of


def flatten(
    table,
    combine: Optional[Dict[str, Callable]] = None,
    split_columns: Optional[Iterable[str]] = None,
):
    """Combine overlapping regions into single rows, similar to bedtools merge."""
    if table.empty:
        return table
    if (table.start.values[1:] >= table.end.cummax().values[:-1]).all():
        return table
    # NB: Input rows and columns should already be sorted like this
    table = table.sort_values(["chromosome", "start", "end"])
    cmb = get_combiners(table, False, combine)
    out = (
        table.groupby(by="chromosome", as_index=False, group_keys=False, sort=False)
        [table.columns]
        .apply(_flatten_overlapping, cmb, split_columns)
        .reset_index(drop=True)
    )
    return out.reindex(
        out.chromosome.apply(sorter_chrom).sort_values(kind="mergesort").index
    )


def _flatten_overlapping(
    table, combine: Dict[str, Callable], split_columns: Optional[Iterable]
):
    """Merge overlapping regions within a chromosome/strand.

    Assume chromosome and (if relevant) strand are already identical, so only
    start and end coordinates are considered.
    """
    if split_columns:
        row_groups = (
            tuple(_flatten_tuples_split(row_group, combine, split_columns))
            for row_group in _nonoverlapping_groups(table, 0)
        )
    else:
        row_groups = (
            tuple(_flatten_tuples(row_group, combine))
            for row_group in _nonoverlapping_groups(table, 0)
        )
    all_row_groups = itertools.chain(*row_groups)
    return pd.DataFrame.from_records(list(all_row_groups), columns=table.columns)


def _flatten_tuples(keyed_rows: Iterable, combine: Dict[str, Callable]):
    """Divide multiple rows where they overlap.

    Parameters
    ----------
    keyed_rows : iterable
        pairs of (non-overlapping-group index, overlapping rows)
    combine : dict
        Mapping of field names to functions applied to combine overlapping
        regions.

    Returns
    -------
    DataFrame
    """
    rows = [kr[1] for kr in keyed_rows]
    first_row = rows[0]
    if len(rows) == 1:
        yield first_row
    else:
        # TODO speed this up! Bottleneck is in dictcomp
        extra_cols = [x for x in first_row._fields[3:] if x in combine]
        breaks = sorted(set(itertools.chain(*[(r.start, r.end) for r in rows])))
        for bp_start, bp_end in zip(breaks[:-1], breaks[1:]):
            # Find the row(s) overlapping this region
            # i.e. any not already seen and not already passed
            rows_in_play = [
                row for row in rows if row.start <= bp_start and row.end >= bp_end
            ]
            # Combine the extra fields of the overlapping regions
            extra_fields = {
                key: combine[key]([getattr(r, key) for r in rows_in_play])
                for key in extra_cols
            }
            yield first_row._replace(start=bp_start, end=bp_end, **extra_fields)


def _flatten_tuples_split(keyed_rows, combine: Dict, split_columns: Optional[Iterable]):
    """Divide multiple rows where they overlap.

    Parameters
    ----------
    keyed_rows : iterable
        pairs of (non-overlapping-group index, overlapping rows)
    combine : dict
        Mapping of field names to functions applied to combine overlapping
        regions.
    split_columns : list or tuple
        Field names where numeric values should be subdivided a region.

    Returns
    -------
    DataFrame
    """
    rows = [kr[1] for kr in keyed_rows]
    first_row = rows[0]
    if len(rows) == 1:
        yield first_row
    else:
        # TODO - use split_columns
        extra_cols = [x for x in first_row._fields[3:] if x in combine]
        breaks = sorted(set(itertools.chain(*[(r.start, r.end) for r in rows])))
        for bp_start, bp_end in zip(breaks[:-1], breaks[1:]):
            # Find the row(s) overlapping this region
            # i.e. any not already seen and not already passed
            rows_in_play = [
                row for row in rows if row.start <= bp_start and row.end >= bp_end
            ]
            # Combine the extra fields of the overlapping regions
            extra_fields = {
                key: combine[key]([getattr(r, key) for r in rows_in_play])
                for key in extra_cols
            }
            yield first_row._replace(start=bp_start, end=bp_end, **extra_fields)


def merge(
    table,
    bp: int = 0,
    stranded: bool = False,
    combine: Optional[Dict[str, Callable]] = None,
):
    """Merge overlapping rows in a DataFrame."""
    if table.empty:
        return table
    gap_sizes = table.start.values[1:] - table.end.cummax().values[:-1]
    if (gap_sizes > -bp).all():
        return table
    if stranded:
        groupkey = ["chromosome", "strand"]
    else:
        # NB: same gene name can appear on alt. contigs
        groupkey = ["chromosome"]
    table = table.sort_values(groupkey + ["start", "end"])
    cmb = get_combiners(table, stranded, combine)
    out = (
        table.groupby(by=groupkey, as_index=False, group_keys=False, sort=False)
        [table.columns]
        .apply(_merge_overlapping, bp, cmb)
        .reset_index(drop=True)
    )
    # Re-sort chromosomes cleverly instead of lexicographically
    return out.reindex(
        out.chromosome.apply(sorter_chrom).sort_values(kind="mergesort").index
    )


def _merge_overlapping(table, bp: int, combine: Dict[str, Callable]):
    """Merge overlapping regions within a chromosome/strand.

    Assume chromosome and (if relevant) strand are already identical, so only
    start and end coordinates are considered.
    """
    merged_rows = [
        _squash_tuples(row_group, combine)
        for row_group in _nonoverlapping_groups(table, bp)
    ]
    return pd.DataFrame.from_records(merged_rows, columns=merged_rows[0]._fields)


def _nonoverlapping_groups(table, bp: int):
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
    group_keys = np.r_[False, gap_sizes > (-bp)].cumsum()
    # NB: pandas groupby seems like the obvious choice over itertools, but it is
    # very slow -- probably because it traverses the whole table (i.e.
    # chromosome) again to select groups, redoing the various inferences that
    # would be worthwhile outside this inner loop. With itertools, we take
    # advantage of the grouping and sorting already done, and don't repeat
    # pandas' traversal and inferences.
    # ENH: Find & use a lower-level, 1-pass pandas function
    keyed_groups = zip(group_keys, table.itertuples(index=False))
    return (row_group for _key, row_group in itertools.groupby(keyed_groups, first_of))


# Squash rows according to a given grouping criterion
# XXX see also segfilter.py
def _squash_tuples(keyed_rows, combine: Dict[str, Callable]):
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
    rows = [kr[1] for kr in keyed_rows]  # list(rows)
    firsttup = rows[0]
    if len(rows) == 1:
        return firsttup
    newfields = {
        key: combiner(pd.Series([getattr(r, key) for r in rows]))
        for key, combiner in combine.items()
    }
    return firsttup._replace(**newfields)
