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
from .combiners import first_of, get_combiners

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

_BF_COLS = ("chromosome", "start", "end")


def flatten(
    table: pd.DataFrame,
    combine: dict[str, Callable] | None = None,
    split_columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Combine overlapping regions into single rows, similar to bedtools merge."""
    if table.empty:
        return table
    if (table.start.to_numpy()[1:] >= table.end.cummax().to_numpy()[:-1]).all():
        return table
    # NB: Input rows and columns should already be sorted like this
    table = table.sort_values(["chromosome", "start", "end"])
    cmb = get_combiners(table, False, combine)
    clustered = bioframe.cluster(
        table, min_dist=0, cols=_BF_COLS, return_cluster_ids=True
    )
    all_rows: list = []
    for _cluster_id, group in clustered.groupby("cluster", sort=False):
        group_rows = group.drop(
            columns=["cluster", "cluster_start", "cluster_end"], errors="ignore"
        )
        all_rows.extend(_flatten_tuples(group_rows, cmb))
    out = pd.DataFrame.from_records(all_rows, columns=table.columns)
    return out.reindex(
        out.chromosome.apply(sorter_chrom).sort_values(kind="mergesort").index
    )


def _flatten_tuples(
    group_df: pd.DataFrame, combine: dict[str, Callable]
) -> list[tuple]:
    """Divide multiple rows where they overlap.

    Split at all coordinate breakpoints and combine extra fields.
    """
    rows = list(group_df.itertuples(index=False))
    if len(rows) == 1:
        return rows
    first_row = rows[0]
    extra_cols = [x for x in first_row._fields[3:] if x in combine]
    breaks = sorted(set(itertools.chain(*[(r.start, r.end) for r in rows])))
    result = []
    for bp_start, bp_end in itertools.pairwise(breaks):
        rows_in_play = [
            row for row in rows if row.start <= bp_start and row.end >= bp_end
        ]
        extra_fields = {
            key: combine[key]([getattr(r, key) for r in rows_in_play])
            for key in extra_cols
        }
        result.append(first_row._replace(start=bp_start, end=bp_end, **extra_fields))
    return result


def merge(
    table: pd.DataFrame,
    bp: int = 0,
    stranded: bool = False,
    combine: dict[str, Callable] | None = None,
) -> pd.DataFrame:
    """Merge overlapping rows in a DataFrame."""
    if table.empty:
        return table
    gap_sizes = table.start.to_numpy()[1:] - table.end.cummax().to_numpy()[:-1]
    if (gap_sizes > -bp).all():
        return table
    if stranded:
        groupkey = ["chromosome", "strand"]
    else:
        groupkey = ["chromosome"]
    table = table.sort_values([*groupkey, "start", "end"])
    cmb = get_combiners(table, stranded, combine)
    min_dist = max(0, -bp)
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
    result_rows = []
    for _cluster_id, group in clustered.groupby("cluster", sort=False):
        if len(group) == 1:
            result_rows.append(group[data_cols].iloc[0])
        else:
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
            result_rows.append(pd.Series(vals))
    out = pd.DataFrame(result_rows, columns=data_cols).reset_index(drop=True)
    # Re-sort chromosomes cleverly instead of lexicographically
    return out.reindex(
        out.chromosome.apply(sorter_chrom).sort_values(kind="mergesort").index
    )


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
    if table.empty or len(table) == 1:
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
        return table

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
    return pd.DataFrame(result_rows, columns=table.columns).reset_index(drop=True)
