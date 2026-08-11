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

Numeric: TypeAlias = int | float | np.number

if TYPE_CHECKING:
    from collections.abc import Callable, Generator, Iterator, Sequence

    from numpy import ndarray

    # One coordinate per interval, whether a table's rows or a query's regions
    Coords: TypeAlias = Sequence[Numeric] | pd.Series | ndarray
    # ... or none at all, for a bound the caller left open
    Bounds: TypeAlias = Coords | None

# Stand-ins for the open side of a one-sided query, chosen to fall outside the
# coordinates rather than at their edge: genomic positions are non-negative,
# and none exceeds what an int64 column can hold.
_BEFORE_FIRST_BASE = -1
_PAST_LAST_BASE = np.iinfo(np.int64).max

# One shared placeholder instead of a fresh array per yield: allocating one
# costs ~200 ns against ~7 us of per-row work, and the missing-chromosome path
# yields one for every row of `other`. Safe to hand out repeatedly because a
# zero-length array has no elements to write.
_EMPTY_SELECTION = np.empty(0, dtype=np.intp)


def point_aware_ends(starts: Coords, ends: Coords) -> ndarray:
    """Give each zero-width interval the one base its start names.

    The overlap predicate, stated once for the whole package: two intervals
    overlap when their half-open coordinate ranges share a base, and an
    interval whose end equals its start is the point at that start -- the BED
    spelling of an insertion site -- rather than the empty set. Rewriting such
    an end as `start + 1` reduces the point to an ordinary one-base interval,
    so every caller needs only the strict half-open comparison.

    This is bioframe's rule: `bioframe.overlap` applies the same rewrite to
    both of its inputs, in `core.arrops._convert_points_to_len1_segments`.
    Reproducing it here is what keeps the `searchsorted` searches below
    answering the same question as the bioframe-backed
    :meth:`GenomicArray.intersection`.

    It governs overlap *selection*; clustering is a different relation, and
    `skgenome.merge` explains why it stays strict.
    """
    start_arr = np.asarray(starts)
    end_arr = np.asarray(ends)
    # Only `end == start` is rewritten. An inverted row, `end < start`, is
    # malformed input that neither half of the package repairs, and
    # `max(end, start + 1)` would quietly rewrite those too -- eight such
    # segments ship in `test/formats/amplicon.cns`. `np.where` evaluates both
    # branches, so reading `start + 1` rather than `end + 1` also keeps
    # `_PAST_LAST_BASE` out of an addition that would wrap in silence.
    return np.where(end_arr == start_arr, start_arr + 1, end_arr)


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
    starts: Bounds,
    ends: Bounds,
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
            # Clip to the region's bounds -- `is not None`, not truthiness. A
            # bound of 0 is a real one: a zero-width query at the origin now
            # selects the rows covering base 0, and a falsy `end_val` left
            # them reaching past it, unclipped. (`clip(lower=0)` is a no-op on
            # non-negative starts, so only the end bound ever did work, but a
            # guard that is right on one side and wrong on the other is worse
            # than either.)
            if start_val is not None:
                subtable.start = subtable.start.clip(lower=start_val)
            if end_val is not None:
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
    starts: Bounds,
    ends: Bounds,
    mode: str,
) -> Generator[tuple[slice | ndarray, Numeric | None, Numeric | None], None, None]:
    """Iterate through sub-ranges.

    `table`'s `start` column must ascend, since that is what lets a binary
    search place the query. Callers that group by chromosome first --
    `by_ranges`, `iter_slices` -- hold that per group, and `iter_ranges` does
    when given a chromosome to select; searching a whole array at once, as
    `in_range(None, ...)` does, needs it end to end. Rows out of order raise
    `ValueError` rather than being answered wrongly.
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

    # Cost: one extra pass over `start`, beside the pass over `end` below that
    # chooses between the two implementations, and the binary searches
    # themselves.
    require_ascending_starts(table)

    n_regions = max(
        0 if starts is None else len(starts), 0 if ends is None else len(ends)
    )
    if not n_regions:
        # No regions to select, so nothing to yield: one region per query, and
        # there are none. (Not the same as the `starts is None and ends is
        # None` case above, which asks for the whole table.)
        return

    # Both implementations below want a start and an end for every region, so
    # stand in for whichever side a one-sided query left open. A real bound of
    # 0 is not the same query as no bound at all -- `(0, 0)` is the point at
    # the origin, while `(None, 0)` reaches only below base 0, where nothing
    # lies -- so the stand-in has to sit outside the coordinates rather than at
    # their edge, where the point rule would widen it into one.
    if starts is None or not len(starts):
        starts = np.full(n_regions, _BEFORE_FIRST_BASE)
    if ends is None or not len(ends):
        ends = np.full(n_regions, _PAST_LAST_BASE)

    # Apply the overlap predicate once, here, so that the two implementations
    # below and the test choosing between them all read the same ends, and so
    # that neither implementation can reach `table.end` itself.
    table_ends = point_aware_ends(table.start, table.end)
    query_starts = np.asarray(starts)
    query_ends = point_aware_ends(query_starts, ends)

    selections: Iterator[slice | ndarray]
    # `pd.Index(...).is_monotonic_increasing` short-circuits, but constructing
    # the Index costs more than the scan it saves: measured 1.6-10x slower
    # than the comparison below over 2.6k-200k rows, sorted or not. Both
    # answer True for an empty or single-row array and False where a
    # coordinate is NaN.
    if bool((table_ends[1:] >= table_ends[:-1]).all()):
        selections = _irange_simple(table, table_ends, query_starts, query_ends, mode)
    else:
        # A bin nested inside another leaves `end` out of ascending order, and
        # `_irange_simple` reads `end` with `searchsorted` too -- to place a
        # start in 'outer' mode, an end in 'inner' mode. The masks
        # `_irange_nested` builds compare `end` elementwise instead.
        selections = _irange_nested(table, table_ends, query_starts, query_ends, mode)
    # Report the bounds as normalised above, not the point-aware widening of
    # them: 'trim' mode clips its output to these, and clipping to `start + 1`
    # would hand back a base the query never asked for. The open-side stand-ins
    # reach it too, where `clip(lower=-1)` and `clip(upper=int64max)` are the
    # no-ops an absent bound should be.
    yield from zip(selections, starts, ends, strict=True)


def require_ascending_starts(table: pd.DataFrame) -> None:
    """Raise `ValueError` unless `table`'s `start` column ascends.

    The precondition every `searchsorted` against a coordinate column carries:
    handed rows in some other order it returns a position unrelated to the
    query, silently, and the caller acts on a row it never asked for. Callers
    that place a coordinate without wanting a row selection back need the
    check on its own, which is why it is not folded into `idx_ranges`.

    `table` needs a `chromosome` column as well as `start`, since that is what
    tells a genuinely unsorted table apart from one merely spanning several
    chromosomes -- which every whole-genome array does, coordinates restarting
    at each one.
    """
    if not table.start.is_monotonic_increasing:
        _raise_unsorted(table)


def _raise_unsorted(table: pd.DataFrame) -> NoReturn:
    """Report what is out of order about `table`, and how to fix it."""
    starts = table["start"]
    chroms = table["chromosome"].unique()
    # Diagnose the cause once: each carries its own remedy, and the widest
    # explanation is the right one to give.
    if len(chroms) > 1:
        # Coordinates restart at every chromosome, so an array holding more
        # than one is out of order by construction. This is a search that
        # named none of them.
        shown = ", ".join(map(str, chroms[:3]))
        opening = "Genomic intervals must be sorted by start position"
        detail = (
            f"these {len(table)} rows span {len(chroms)} chromosomes "
            f"({shown}{', ...' if len(chroms) > 3 else ''})"
        )
        remedy = "name the chromosome to search rather than passing None"
    elif starts.isna().any():
        # A missing position compares false against every other, so the column
        # counts as unordered with no row behind the one before it.
        opening = "Genomic intervals must have a start position"
        detail = f"start is missing in {int(starts.isna().sum())} of {len(table)} rows"
        remedy = "drop the rows whose start is missing"
    else:
        # Nothing is missing, so some row's start is simply behind the one
        # before it. Compare the column against itself shifted rather than
        # differencing it, so that reporting the problem cannot fail on
        # whatever the caller put in the column.
        values = starts.to_numpy()
        row = int(np.flatnonzero(values[1:] < values[:-1])[0]) + 1
        opening = "Genomic intervals must be sorted by start position"
        detail = (
            f"start={values[row]} follows start={values[row - 1]}, at row "
            f"{row} of {len(table)}"
        )
        remedy = "sort the rows by start -- GenomicArray.sort, or sort_values"
    raise ValueError(f"{opening} to be searched: {detail}. To fix, {remedy}.")


def _irange_simple(
    table: pd.DataFrame,
    table_ends: ndarray,
    starts: ndarray,
    ends: ndarray,
    mode: str,
) -> Iterator[slice]:
    """Select each region by slice, for a table whose bins do not nest.

    Both coordinate columns ascend, so each region's rows are one run of the
    table and binary search finds its two ends. `table_ends` and `ends` carry
    the overlap predicate, `point_aware_ends`; `table.end` is not read here.
    """
    if mode == "inner":
        # Only the rows the region encloses
        start_idxs = table.start.searchsorted(starts)
        end_idxs = table_ends.searchsorted(ends, "right")
    else:
        # Also the rows straddling either boundary
        start_idxs = table_ends.searchsorted(starts, "right")
        end_idxs = table.start.searchsorted(ends)

    for start_idx, end_idx in zip(start_idxs, end_idxs, strict=True):
        yield slice(start_idx, end_idx)


def _irange_nested(
    table: pd.DataFrame,
    table_ends: ndarray,
    starts: ndarray,
    ends: ndarray,
    mode: str,
) -> Iterator[ndarray]:
    """Select each region by mask, for a table whose bins nest.

    A nested bin leaves `end` out of ascending order, so only `start` can be
    searched; `end` is compared row by row instead. `table_ends` and `ends`
    carry the overlap predicate, `point_aware_ends`; `table.end` is not read
    here.
    """
    # ENH: Binary Interval Search (BITS) or Layer&Quinlan(2015)
    for start_val, end_val in zip(starts, ends, strict=True):
        # Mask of table rows to keep for this query region
        if mode == "inner":
            # Only the rows the region encloses
            region_mask = table_ends <= end_val
            region_mask[: int(table.start.searchsorted(start_val))] = 0
        else:
            # Also the rows straddling either boundary
            region_mask = table_ends > start_val
            region_mask[int(table.start.searchsorted(end_val)) :] = 0

        yield region_mask
