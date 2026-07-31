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
from pandas.api.types import is_integer_dtype

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
        total is preserved instead of being duplicated into every piece. An
        integer column stays integral, rounded piece by piece, since a count
        such as ``probes`` is not meaningful as a fraction. Pass an empty
        sequence to combine these columns like any other.

        Conservation holds for a column with an additive combiner, the default
        for both of these, and for regions of positive length: a zero-width
        row contributes a breakpoint but has no sub-interval to carry its
        value, so it is dropped as it always has been.

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
    # Filtered to the columns that have a combiner, since a column with none is
    # taken from its covering row verbatim and so cannot be apportioned
    split_cols = {
        col
        for col in (_SPLIT_COLUMNS if split_columns is None else split_columns)
        if col in cmb
    }
    # NB: Input rows and columns should already be sorted like this
    table = table.sort_values(["chromosome", "start", "end"])
    # Bookended rows share no breakpoint to split at, so cluster only the rows
    # that genuinely overlap -- the same distinction `_nothing_to_cluster`
    # makes above when it lets a whole table through.
    clustered = bioframe.cluster(
        table, min_dist=None, cols=_BF_COLS, return_cluster_ids=True
    )
    data_cols = _data_columns(clustered)
    out = _splice_clusters(
        clustered,
        data_cols,
        lambda group: _flatten_cluster(group[data_cols], cmb, split_cols),
    )
    # A count stays a count: 'probes' is documented as a number of bins, and
    # `export vcf` reads a fractional probe count as a corrupt segment and drops
    # it. Rounding the shares back costs under one count per piece.
    for col in split_cols:
        if is_integer_dtype(table[col]):
            out[col] = out[col].round().astype(table[col].dtype)
    # Re-sort chromosomes cleverly instead of lexicographically
    out = out.reindex(
        out.chromosome.apply(sorter_chrom).sort_values(kind="mergesort").index
    )
    return _fill_unnamed(out, cmb)


def _data_columns(clustered: pd.DataFrame) -> list[str]:
    """The columns of a `bioframe.cluster` result that came from the input."""
    return [
        col
        for col in clustered.columns
        if col not in ("cluster", "cluster_start", "cluster_end")
    ]


def _splice_clusters(
    clustered: pd.DataFrame,
    data_cols: list[str],
    expand: Callable[[pd.DataFrame], list[dict]],
) -> pd.DataFrame:
    """Rebuild a clustered table, replacing each multi-row cluster with `expand`.

    Rows that cluster alone are emitted verbatim; only the rows that really
    cluster together are worth a Python-level pass. When a bait file holds a
    handful of overlaps nearly every cluster is a singleton, and visiting them
    all otherwise dominates the runtime.
    """
    cluster_ids = clustered["cluster"]
    is_first = ~cluster_ids.duplicated().to_numpy()
    in_multi_row_cluster = cluster_ids.duplicated(keep=False).to_numpy()
    out = clustered.loc[is_first, data_cols].reset_index(drop=True)
    if not in_multi_row_cluster.any():
        return out
    # `groupby(sort=False)` yields clusters in order of first appearance, which
    # is the order `is_first` picked them out of `out`. The rows replacing a
    # cluster take the index of the one row standing in for it there, and
    # `mergesort` is stable, so they keep the order `expand` gave them.
    changing = in_multi_row_cluster[is_first]
    rows: list[dict] = []
    index: list[int] = []
    for position, (_cluster_id, group) in zip(
        np.flatnonzero(changing).tolist(),
        clustered[in_multi_row_cluster].groupby("cluster", sort=False),
        strict=True,
    ):
        new_rows = expand(group)
        rows.extend(new_rows)
        index.extend(itertools.repeat(position, len(new_rows)))
    # Splice via `concat` so pandas widens each column's dtype to fit the
    # combined values, e.g. a joined name landing in an all-NaN float column.
    replacement = pd.DataFrame(rows, columns=data_cols, index=index)
    return (
        pd.concat([out[~changing], replacement])
        .sort_index(kind="mergesort")
        .reset_index(drop=True)
    )


def _flatten_cluster(
    group_df: pd.DataFrame, cmb: dict[str, Callable], split_cols: set[str]
) -> list[dict]:
    """Divide one cluster of overlapping rows at each of their breakpoints.

    Every sub-interval draws its values from the rows covering that
    sub-interval alone: the columns in `cmb` through their combiner, the rest
    from the first covering row. Columns in `split_cols` are apportioned by
    length before being combined.
    """
    # Address the rows positionally: `itertuples` renames any column whose name
    # is not an identifier, and a tabular input's header is the user's own
    positions = {col: i for i, col in enumerate(group_df.columns)}
    i_start = positions["start"]
    i_end = positions["end"]
    rows = list(group_df.itertuples(index=False, name=None))
    breaks = sorted({bound for row in rows for bound in (row[i_start], row[i_end])})
    result = []
    for bp_start, bp_end in itertools.pairwise(breaks):
        # The cluster spans a contiguous range and no breakpoint falls strictly
        # within this sub-interval, so at least one row covers all of it
        covering = [
            row for row in rows if row[i_start] <= bp_start and row[i_end] >= bp_end
        ]
        piece = {}
        for col, i_col in positions.items():
            if col == "start":
                piece[col] = bp_start
            elif col == "end":
                piece[col] = bp_end
            elif col in cmb:
                if col in split_cols:
                    # A row's share of the quantity is its share of the length
                    piece[col] = cmb[col](
                        [
                            row[i_col]
                            * (bp_end - bp_start)
                            / (row[i_end] - row[i_start])
                            for row in covering
                        ]
                    )
                else:
                    piece[col] = cmb[col]([row[i_col] for row in covering])
            else:
                piece[col] = covering[0][i_col]
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
    data_cols = _data_columns(clustered)

    def combine_cluster(group: pd.DataFrame) -> list[dict]:
        """Combine one cluster's rows into the single row of their union."""
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
        return [vals]

    out = _splice_clusters(clustered, data_cols, combine_cluster)
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
