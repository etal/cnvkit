"""DataFrame-level merging operations.

Merge overlapping regions into single rows, similar to bedtools merge.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import bioframe
import numpy as np
import pandas as pd

from .chromsort import sorter_chrom
from .combiners import first_of, get_combiners, join_strings

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

_BF_COLS = ("chromosome", "start", "end")
# Columns that `flatten` apportions across the pieces of a split region, being
# quantities spread over its length rather than properties of it
_SPLIT_COLUMNS = ("weight", "probes")


def _fill_unnamed(table: pd.DataFrame, cmb: dict[str, Callable]) -> pd.DataFrame:
    """Fill non-string values in string-join columns with the "-" placeholder.

    ``merge``/``flatten``/``squash`` emit singleton (non-overlapping) rows
    verbatim, so a NaN-float gene/accession from a ragged input survives
    unjoined. Bring those into line with the multi-row output of
    :func:`join_strings`, which renders missing names as ``"-"``. Only columns
    that actually use ``join_strings`` are touched, so caller-supplied combiners
    are respected. A no-op returning the input unchanged when every value is
    already a string, so well-formed data is neither copied nor altered.
    """
    fix = {}
    for col, fn in cmb.items():
        if fn is join_strings and col in table.columns:
            nonstr = ~table[col].map(lambda v: isinstance(v, str)).to_numpy()
            if nonstr.any():
                fix[col] = nonstr
    if not fix:
        return table
    table = table.copy()
    for col, nonstr in fix.items():
        table.loc[nonstr, col] = "-"
    return table


def _nothing_to_cluster(table: pd.DataFrame, bp: int) -> bool:
    """Whether no two rows would cluster together, so the table passes through.

    Rows cluster when the gap between them is at most ``-bp`` bases: ``bp = 0``
    clusters bookended rows, and ``bp = 1`` clusters only rows that genuinely
    overlap.

    Each row is compared against the running maximum end of the preceding rows
    on the same chromosome, which assumes the table is sorted -- as everything
    read through :mod:`skgenome.tabio` is. Tables whose rows are out of order,
    or whose chromosomes are not in contiguous blocks, fall through to the
    clustering path, which sorts.
    """
    chromosomes = table.chromosome.to_numpy()
    changes = chromosomes[1:] != chromosomes[:-1]
    if 1 + int(changes.sum()) != table.chromosome.nunique():
        # A chromosome occurs in more than one block, so the running maximum
        # below would skip comparisons between rows that are not adjacent
        return False
    covered = table.groupby("chromosome", sort=False).end.cummax().to_numpy()
    gap_sizes = table.start.to_numpy()[1:] - covered[:-1]
    return bool((changes | (gap_sizes > -bp)).all())


def flatten(
    table: pd.DataFrame,
    combine: dict[str, Callable] | None = None,
    split_columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Split overlapping regions at every breakpoint, combining extra fields.

    Unlike :func:`merge`, which emits an overlapping group's union as one row,
    this emits one row per distinct sub-interval, so partially overlapping
    regions yield fragments narrower than either input.

    Parameters
    ----------
    table : DataFrame
        Genomic intervals, with at least chromosome, start and end columns.
        Must be sorted by chromosome, start, end.
    combine : dict, optional
        Column-to-function mapping for aggregating the fields of the rows that
        cover each sub-interval. See :func:`skgenome.combiners.get_combiners`.
    split_columns : iterable of str, optional
        Columns holding a quantity spread over the length of a region rather
        than a property of it -- ``weight`` and ``probes`` by default. Each
        row's value is apportioned among its sub-intervals in proportion to
        their lengths before the covering rows are combined, so the column's
        total is preserved instead of being duplicated into every piece. Pass
        an empty sequence to combine these columns like any other; columns
        without a combiner are ignored.

    Returns
    -------
    DataFrame
        A new table with the overlapping rows split at each breakpoint.
    """
    if table.empty:
        return table
    cmb = get_combiners(table, False, combine)
    if _nothing_to_cluster(table, 1):
        return _fill_unnamed(table, cmb)
    split_cols = {
        col
        for col in (_SPLIT_COLUMNS if split_columns is None else split_columns)
        if col in cmb
    }
    # NB: Input rows and columns should already be sorted like this
    table = table.sort_values(["chromosome", "start", "end"])
    clustered = bioframe.cluster(
        table, min_dist=0, cols=_BF_COLS, return_cluster_ids=True
    )
    all_rows: list = []
    for _cluster_id, group in clustered.groupby("cluster", sort=False):
        group_rows = group.drop(
            columns=["cluster", "cluster_start", "cluster_end"], errors="ignore"
        )
        all_rows.extend(_flatten_cluster(group_rows, cmb, split_cols))
    out = pd.DataFrame(all_rows, columns=table.columns)
    out = out.reindex(
        out.chromosome.apply(sorter_chrom).sort_values(kind="mergesort").index
    )
    return _fill_unnamed(out, cmb)


def _flatten_cluster(
    group_df: pd.DataFrame, combine: dict[str, Callable], split_cols: set[str]
) -> list[dict]:
    """Divide one cluster of overlapping rows at each of their breakpoints.

    Every sub-interval draws its values from the rows covering that
    sub-interval alone: the columns in `combine` through their combiner, the
    rest from the first covering row. Columns in `split_cols` are apportioned
    by length before being combined.
    """
    columns = list(group_df.columns)
    rows = list(group_df.itertuples(index=False))
    breaks = sorted({bound for row in rows for bound in (row.start, row.end)})
    result = []
    for bp_start, bp_end in itertools.pairwise(breaks):
        # The cluster spans a contiguous range and no breakpoint falls strictly
        # within this sub-interval, so at least one row covers all of it
        covering = [row for row in rows if row.start <= bp_start and row.end >= bp_end]
        piece = {}
        for col in columns:
            if col == "start":
                piece[col] = bp_start
            elif col == "end":
                piece[col] = bp_end
            elif col in combine:
                values = [getattr(row, col) for row in covering]
                if col in split_cols:
                    values = [
                        val * (bp_end - bp_start) / (row.end - row.start)
                        for val, row in zip(values, covering, strict=True)
                    ]
                piece[col] = combine[col](values)
            else:
                piece[col] = getattr(covering[0], col)
        result.append(piece)
    return result


def merge(
    table: pd.DataFrame,
    bp: int = 0,
    stranded: bool = False,
    combine: dict[str, Callable] | None = None,
) -> pd.DataFrame:
    """Merge rows that overlap, or that lie within a given gap of each other.

    Parameters
    ----------
    table : DataFrame
        Genomic intervals, with at least chromosome, start and end columns.
        Must be sorted by chromosome, start, end.
    bp : int
        Rows are merged when the gap between them is at most ``-bp`` bases.
        The default 0 merges bookended rows too, as ``bedtools merge`` does,
        and negative values bridge gaps of that many bases. Any value above 0
        merges only rows that genuinely overlap, leaving bookended rows -- the
        consecutive tiles of a bait file, say -- separate; bioframe cannot
        express a finer distinction, so 1 is the effective maximum.
    stranded : bool
        Merge only rows on the same strand, keeping the strand in the output.
    combine : dict, optional
        Column-to-function mapping for aggregating the merged rows' other
        fields. See :func:`skgenome.combiners.get_combiners`.

    Returns
    -------
    DataFrame
        A new table with the overlapping (or near-enough) rows merged.
    """
    if table.empty:
        return table
    cmb = get_combiners(table, stranded, combine)
    # bioframe clusters rows at most `min_dist` bases apart, and with None only
    # rows that truly overlap -- the finest distinction it can express, so any
    # `bp` above 1 is the same as 1.
    bp = min(bp, 1)
    if _nothing_to_cluster(table, bp):
        return _fill_unnamed(table, cmb)
    groupkey = ["chromosome", "strand"] if stranded else ["chromosome"]
    table = table.sort_values([*groupkey, "start", "end"])
    min_dist = None if bp == 1 else -bp
    on = ["strand"] if stranded else None
    clustered = bioframe.cluster(
        table,
        min_dist=min_dist,
        cols=_BF_COLS,
        on=on,
        return_cluster_ids=True,
        return_cluster_intervals=True,
    )
    data_cols = [
        c
        for c in clustered.columns
        if c not in ("cluster", "cluster_start", "cluster_end")
    ]
    # Rows that cluster alone are emitted verbatim; only the rows that really
    # combine are worth a Python-level pass. When a bait file holds a handful
    # of duplicates nearly every cluster is a singleton, and visiting them all
    # otherwise dominates the runtime.
    cluster_ids = clustered["cluster"]
    is_first = ~cluster_ids.duplicated().to_numpy()
    in_multi_row_cluster = cluster_ids.duplicated(keep=False).to_numpy()
    out = clustered.loc[is_first, data_cols].reset_index(drop=True)
    if in_multi_row_cluster.any():
        combined_rows = []
        for _cluster_id, group in clustered[in_multi_row_cluster].groupby(
            "cluster", sort=False
        ):
            vals = {}
            for col in data_cols:
                if col == "start":
                    vals[col] = group["cluster_start"].iloc[0]
                elif col == "end":
                    vals[col] = group["cluster_end"].iloc[0]
                elif col in cmb:
                    vals[col] = cmb[col](group[col].values)
                else:
                    vals[col] = group[col].iloc[0]
            combined_rows.append(vals)
        # `groupby(sort=False)` yields clusters in order of first appearance,
        # which is the order `is_first` picked them out of `out`. Splice via
        # `concat` so pandas widens each column's dtype to fit the combined
        # values, e.g. a joined name landing in an all-NaN float column.
        combining = in_multi_row_cluster[is_first]
        replacement = pd.DataFrame(
            combined_rows, columns=data_cols, index=np.flatnonzero(combining)
        )
        out = pd.concat([out[~combining], replacement]).sort_index(kind="mergesort")
    # Re-sort chromosomes cleverly instead of lexicographically
    out = out.reindex(
        out.chromosome.apply(sorter_chrom).sort_values(kind="mergesort").index
    )
    return _fill_unnamed(out, cmb)


def squash(
    table: pd.DataFrame,
    by: str | None = None,
    combine: dict[str, Callable] | None = None,
) -> pd.DataFrame:
    """Combine consecutive adjacent rows into single rows.

    Parameters
    ----------
    table : DataFrame
        Genomic intervals with at least chromosome, start, end columns.
        Must be sorted by chromosome, start, end.
    by : str or None
        If given, only combine consecutive rows that have the same value
        in this column (e.g. ``"gene"``).  If None, combine all adjacent
        rows on the same chromosome.
    combine : dict or None
        Column-to-function mappings for aggregation.  See
        :func:`get_combiners`.

    Returns
    -------
    DataFrame
        A new DataFrame with consecutive adjacent rows combined.
    """
    if table.empty:
        return table

    cmb = get_combiners(table, stranded=False, combine=combine)
    if by is not None:
        # Values are identical within each group, so just take the first
        cmb[by] = first_of

    # Vectorised adjacency detection
    chroms = table["chromosome"].to_numpy()
    starts = table["start"].to_numpy()
    ends = table["end"].to_numpy()
    is_adjacent = (chroms[:-1] == chroms[1:]) & (ends[:-1] == starts[1:])
    if by is not None:
        by_vals = table[by].to_numpy()
        is_adjacent &= by_vals[:-1] == by_vals[1:]

    # Nothing to squash -- short-circuit
    if not is_adjacent.any():
        return _fill_unnamed(table, cmb)

    # Assign group IDs: increment when rows are *not* adjacent
    group_ids = np.concatenate([[0], np.cumsum(~is_adjacent)])

    result_rows: list[pd.Series] = []
    for _gid, group in table.groupby(group_ids, sort=False):
        if len(group) == 1:
            result_rows.append(group.iloc[0])
        else:
            vals: dict = {}
            for col in table.columns:
                if col == "start":
                    vals[col] = group["start"].iloc[0]
                elif col == "end":
                    vals[col] = group["end"].iloc[-1]
                elif col in cmb:
                    vals[col] = cmb[col](group[col].values)
                else:
                    vals[col] = group[col].iloc[0]
            result_rows.append(pd.Series(vals))
    out = pd.DataFrame(result_rows, columns=table.columns).reset_index(drop=True)
    return _fill_unnamed(out, cmb)
