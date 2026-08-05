"""DataFrame-level intersection operations.

Calculate overlapping regions, similar to bedtools intersect.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""

from __future__ import annotations

from itertools import repeat
from typing import TYPE_CHECKING, Any, NoReturn, TypeAlias

import numpy as np
import pandas as pd

from .combiners import first_of, join_strings, make_const

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterator, Sequence

    from numpy import ndarray

Numeric: TypeAlias = int | float | np.number

# One shared placeholder instead of a fresh array per yield: allocating one
# costs ~200 ns against ~7 us of per-row work, and the missing-chromosome path
# yields one for every row of `other`. Safe to hand out repeatedly because a
# zero-length array has no elements to write.
_EMPTY_SELECTION = np.empty(0, dtype=np.intp)


def by_ranges(
    table: pd.DataFrame, other: pd.DataFrame, mode: str, keep_empty: bool
) -> Generator[tuple, None, None]:
    """Group rows by another GenomicArray's bin coordinate ranges.

    Yields one (row of `other`, overlapping rows of `table`) pair per row of
    `other`, in `other`'s row order, unless `keep_empty` is false and the
    selection is empty.
    """
    # Build one sub-iterator per chromosome and draw from them in `other`'s row
    # order; see `iter_slices` for why the grouped order of `by_shared_chroms`
    # will not do.
    per_chrom: dict[Any, Iterator] = {}
    for chrom, bin_rows, src_rows in by_shared_chroms(other, table):
        if src_rows is None or not len(src_rows):
            # `table` has no rows on this chromosome, so every selection is
            # empty. `by_shared_chroms` yields no empty group today; the length
            # check is defensive, since `idx_ranges` would collapse an empty
            # table to a single slice(None) and silently under-yield.
            # ENH: empty dframe matching table, rather than a bare sequence
            per_chrom[chrom] = repeat((), len(bin_rows))
        else:
            per_chrom[chrom] = iter_ranges(
                src_rows, None, bin_rows["start"], bin_rows["end"], mode
            )
    for bin_row in other.itertuples(index=False):
        subrange = next(per_chrom[bin_row.chromosome])
        if keep_empty or len(subrange):
            yield bin_row, subrange


def by_shared_chroms(
    table: pd.DataFrame, other: pd.DataFrame
) -> Iterator[tuple[str, pd.DataFrame, pd.DataFrame | None]]:
    """Group rows for both `table` and `other` by matching chromosome names.

    Yields one group per chromosome of `table`, in order of first appearance;
    the third element is None where `other` has no rows on that chromosome.
    """
    # When both `table` and `other` contain only one chromosome each, and it's
    # the same chromosome, we can just return the original tables.
    table_chr, other_chr = set(table["chromosome"]), set(other["chromosome"])
    if len(table_chr) == 1 and table_chr == other_chr:
        yield table["chromosome"].iat[0], table, other
    else:
        # C416 would suggest `dict(...)`, but `dict()` on a pandas
        # DataFrameGroupBy follows the mapping protocol (via __getitem__),
        # not the iterable-of-pairs protocol, so it fails at runtime with
        # `'str' object is not callable`. Keep the comprehension.
        other_chroms = {c: o for c, o in other.groupby("chromosome", sort=False)}  # noqa: C416
        for chrom, ctable in table.groupby("chromosome", sort=False):
            yield chrom, ctable, other_chroms.get(chrom)


def into_ranges(
    source: pd.DataFrame,
    dest: pd.DataFrame,
    src_col: str,
    default: Any,
    summary_func: Callable | None,
) -> pd.Series:
    """Group a column in `source` by regions in `dest` and summarize."""
    if not len(source) or not len(dest):
        return pd.Series(np.repeat(default, len(dest)), index=dest.index)

    if summary_func is None:
        # Choose a type-appropriate summary function
        elem = source[src_col].iat[0]
        if isinstance(elem, str | np.bytes_):
            summary_func = join_strings
        elif isinstance(elem, float | np.float64):
            summary_func = np.nanmedian
        else:
            summary_func = first_of
    elif not callable(summary_func):  # type: ignore[unreachable]
        # Just fill in the given value, I suppose.
        summary_func = make_const(summary_func)  # type: ignore[unreachable]

    def series2value(ser):
        if len(ser) == 0:
            return default
        if len(ser) == 1:
            return ser.iat[0]
        return summary_func(ser)

    column = source[src_col]
    result = [
        series2value(column[slc]) for slc in iter_slices(source, dest, "outer", True)
    ]
    # Index by `dest`'s labels: callers assign the result back onto `dest`
    # columns, and pandas aligns that assignment by label, not by position.
    return pd.Series(result, index=dest.index)


def iter_ranges(
    table: pd.DataFrame,
    chrom: str | None,
    starts: Sequence[Numeric] | None,
    ends: Sequence[Numeric] | None,
    mode: str,
) -> Iterator[pd.DataFrame]:
    """Iterate through sub-ranges."""
    assert mode in ("inner", "outer", "trim")
    # Optional if we've already subsetted by chromosome (not checked!)
    if chrom:
        assert isinstance(chrom, str)  # ENH: accept array?
        try:
            table = table[table["chromosome"] == chrom]
        except KeyError as exc:
            raise KeyError(f"Chromosome {chrom} is not in this probe set") from exc
    for region_idx, start_val, end_val in idx_ranges(
        table, starts, ends, "inner" if mode == "inner" else "outer"
    ):
        subtable = table.iloc[region_idx]
        if mode == "trim":
            subtable = subtable.copy()
            # Update 5' endpoints to the boundary
            if start_val:
                subtable.start = subtable.start.clip(lower=start_val)
            # Update 3' endpoints to the boundary
            if end_val:
                subtable.end = subtable.end.clip(upper=end_val)
        yield subtable


def iter_slices(
    table: pd.DataFrame, other: pd.DataFrame, mode: str, keep_empty: bool
) -> Iterator[ndarray]:
    """Yields indices to extract ranges from `table`, in `other`'s row order.

    Returns an iterable of integer arrays that can apply to Series objects,
    i.e. columns of `table`. These indices are of the DataFrame/Series' Index,
    not array coordinates -- so be sure to use DataFrame.loc, Series.loc, or
    Series getitem, as opposed to .iloc or indexing directly into Numpy arrays.

    The yields are always in `other`'s row order. With `keep_empty` there is
    additionally exactly one per row, so callers may pair them positionally;
    without it, a row whose selection is empty is skipped.
    """
    # `by_shared_chroms` groups `other` by chromosome, so consuming it directly
    # would yield in chromosome-of-first-appearance order, which differs from
    # row order once a chromosome's rows are non-contiguous; a caller pairing
    # the yields with rows positionally would then write each row's values onto
    # another chromosome's row. Instead, build one lazy sub-iterator per
    # chromosome and draw from them in `other`'s row order.
    per_chrom: dict[Any, Iterator[ndarray]] = {}
    for chrom, bin_rows, src_rows in by_shared_chroms(other, table):
        if src_rows is None or not len(src_rows):
            # `table` has no rows on this chromosome, so every selection is
            # empty. `by_shared_chroms` yields no empty group today; the length
            # check is defensive, since `idx_ranges` would collapse an empty
            # table to a single slice(None) and silently under-yield.
            per_chrom[chrom] = repeat(_EMPTY_SELECTION, len(bin_rows))
        else:
            per_chrom[chrom] = _chrom_slices(src_rows, bin_rows, mode)
    for chrom in other["chromosome"].to_numpy():
        indices = next(per_chrom[chrom])
        if keep_empty or len(indices):
            yield indices


def _chrom_slices(
    src_rows: pd.DataFrame, bin_rows: pd.DataFrame, mode: str
) -> Iterator[ndarray]:
    """Yield `src_rows` index labels selected by each row of `bin_rows`.

    A named function rather than an inline generator expression: a genexp
    evaluates only its outermost iterable eagerly, so its body would re-read
    the caller's loop variables and every chromosome would slice the last one's
    rows.
    """
    for slc, _start, _end in idx_ranges(src_rows, bin_rows.start, bin_rows.end, mode):
        yield src_rows.index[slc].to_numpy()


def idx_ranges(
    table: pd.DataFrame,
    starts: list[int] | pd.Series,
    ends: list[int] | pd.Series,
    mode: str,
) -> Generator[tuple, None, None]:
    """Iterate through sub-ranges.

    `table` must hold one chromosome's rows in ascending `start` order, as
    every caller here arranges: `by_ranges` and `iter_slices` group by
    chromosome, `iter_ranges` selects one, and `tabio.read` sorts what it
    loads. Rows out of that order raise `ValueError`.
    """
    assert mode in ("inner", "outer")
    # Edge cases: when the `table` is either empty, or both `starts` and `ends`
    # are None, we want to signal the calling function to use the entire table.
    # To do this, we return slice(None), which, when passed to either .loc or
    # .iloc, will do just this. We cannot pass table.index to accomplish this
    # because it will not work with .iloc if the table is already subset by
    # chromosome.
    if not len(table) or (starts is None and ends is None):
        yield slice(None), None, None
        return

    # Both implementations below place their slice boundaries with
    # `searchsorted`, which reads its column as ascending: handed rows in some
    # other order it returns a position unrelated to the query, silently, and
    # the caller gets bins that do not overlap the region it asked about. One
    # pass to rule that out, against the pass over `end` just below and the
    # binary searches themselves.
    if not table.start.is_monotonic_increasing:
        _raise_unsorted(table)

    n_regions = max(
        0 if starts is None else len(starts), 0 if ends is None else len(ends)
    )
    if not n_regions:
        # No regions to select, so nothing to yield: one region per query, and
        # there are none. (Not the same as the `starts is None and ends is
        # None` case above, which asks for the whole table.)
        return

    if table.end.is_monotonic_increasing:
        yield from _irange_simple(table, starts, ends, mode)
        return

    # At least one bin is fully nested in another, leaving `end` out of
    # ascending order, and `_irange_simple` reads `end` with `searchsorted`
    # too -- to place a start in 'outer' mode, an end in 'inner' mode. The
    # masks `_irange_nested` builds compare `end` elementwise instead, so they
    # do not care, but it wants a start and an end for every region: fill in
    # whichever side a one-sided query left open.
    if starts is None or not len(starts):
        starts = np.zeros(n_regions, dtype=np.int_)
    if ends is None or not len(ends):
        ends = [None] * n_regions
    yield from _irange_nested(table, starts, ends, mode)


def _raise_unsorted(table: pd.DataFrame) -> NoReturn:
    """Report the first row out of order, and the likely reason for it."""
    row = int(np.flatnonzero(np.diff(table["start"].to_numpy()) < 0)[0]) + 1
    chroms = table["chromosome"].unique()
    if len(chroms) > 1:
        names = ", ".join(map(str, chroms[:3]))
        remedy = (
            f"this selection spans {len(chroms)} chromosomes "
            f"({names}{', ...' if len(chroms) > 3 else ''}), so name the one "
            "to search rather than passing None"
        )
    else:
        remedy = "sort the array with .sort() before selecting from it"
    raise ValueError(
        "Genomic intervals must be sorted by start position to be searched: "
        f"start={table['start'].iat[row]} follows "
        f"start={table['start'].iat[row - 1]}, at row {row} of {len(table)}. "
        f"To fix, {remedy}."
    )


def _irange_simple(
    table: pd.DataFrame, starts: pd.Series, ends: pd.Series, mode: str
) -> Iterator[tuple[slice, int, int]]:
    """Slice subsets of table when regions are not nested."""
    if starts is not None and len(starts):
        if mode == "inner":
            # Only rows entirely after the start point
            start_idxs = table.start.searchsorted(starts)
        else:
            # Include all rows overlapping the start point
            start_idxs = table.end.searchsorted(starts, "right")
    else:
        # `idx_ranges` has ruled out both sides being absent, so `ends` is
        # there to take the region count from.
        starts = np.zeros(len(ends), dtype=np.int_)
        start_idxs = starts.copy()

    if ends is not None and len(ends):
        if mode == "inner":
            end_idxs = table.end.searchsorted(ends, "right")
        else:
            end_idxs = table.start.searchsorted(ends)
    else:
        end_idxs = np.repeat(len(table), len(starts))
        ends = [None] * len(starts)

    for start_idx, start_val, end_idx, end_val in zip(
        start_idxs, starts, end_idxs, ends, strict=True
    ):
        yield (slice(start_idx, end_idx), start_val, end_val)


def _irange_nested(
    table: pd.DataFrame,
    starts: Sequence | ndarray,
    ends: Sequence | ndarray,
    mode: str,
) -> Iterator[tuple[ndarray, int, int]]:
    """Slice subsets of table when regions are nested."""
    # ENH: Binary Interval Search (BITS) or Layer&Quinlan(2015)
    assert len(starts) == len(ends) > 0
    for start_val, end_val in zip(starts, ends, strict=True):
        # Mask of table rows to keep for this query region
        region_mask = np.ones(len(table), dtype=np.bool_)
        if start_val:
            if mode == "inner":
                # Only rows entirely after the start point
                start_idx = table.start.searchsorted(start_val)
                region_mask[: int(start_idx)] = 0
            else:
                # Include all rows overlapping the start point
                region_mask = table.end.to_numpy() > start_val
        if end_val is not None:
            if mode == "inner":
                # Only rows up to the end point
                region_mask &= table.end.to_numpy() <= end_val
            else:
                # Include all rows overlapping the end point
                end_idx = table.start.searchsorted(end_val)
                region_mask[int(end_idx) :] = 0

        yield region_mask, start_val, end_val
