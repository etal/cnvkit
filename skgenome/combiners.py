"""Combiner functions for Python list-like input."""

from collections.abc import Callable, Iterable, Sequence
from typing import Any

import pandas as pd


def get_combiners(
    table: pd.DataFrame,
    stranded: bool = False,
    combine: dict[str, Callable] | None = None,
) -> dict[str, Callable]:
    """Get a `combine` lookup suitable for `table`.

    Parameters
    ----------
    table : DataFrame
    stranded : bool
    combine : dict or None
        Column names to their value-combining functions, replacing or in
        addition to the defaults.

    Returns
    -------
    dict:
        Column names to their value-combining functions.
    """
    cmb = {
        "chromosome": first_of,
        "start": first_of,
        "end": max,
        "gene": join_strings,
        "accession": join_strings,
        "weight": sum,
        "probes": sum,
    }
    if combine:
        cmb |= combine
    if "strand" not in cmb:
        cmb["strand"] = first_of if stranded else merge_strands
    return {k: v for k, v in cmb.items() if k in table.columns}  # type: ignore[misc]


def first_of(elems: Sequence) -> Any:
    """Return the first element of the input."""
    return elems[0]


def last_of(elems: Sequence) -> Any:
    """Return the last element of the input."""
    return elems[-1]


max_of = max


def join_strings(
    elems: Iterable,
    sep: str = ",",
    ignore: tuple[str, ...] = (),
) -> str:
    """Join a Series of unique names, skipping NaN values and ignored names.

    Inputs may themselves already be joined labels -- a bin spanning two genes
    carries the label ``"GENEA,GENEB"`` -- and combining operations chain:
    :func:`skgenome.merge.merge` joins the labels of overlapping baits, and
    :meth:`GenomicArray.into_ranges` then joins those already-joined labels
    again when re-labelling the bins cut from them. Names are therefore
    deduplicated across the separator, so joining is idempotent and no name is
    repeated in the result.

    Parameters
    ----------
    elems : iterable
        Values to join. Non-string elements (e.g. NaN) are silently skipped.
    sep : str
        Separator between joined names, and within any already-joined input.
    ignore : tuple of str
        Names to exclude from the result (e.g. placeholder gene names).

    Returns
    -------
    str
        The joined string, or ``"-"`` if no valid names remain.
    """
    unique_strs = dict.fromkeys(
        name
        for e in elems
        if isinstance(e, str)
        for name in e.split(sep)
        if name not in ignore
    )
    return sep.join(unique_strs) or "-"


def merge_strands(elems: Sequence) -> str:
    """Summarize the given strands as '+', '-', or '.' (both/mixed)"""
    strands = set(elems)
    if len(strands) > 1:
        return "."
    return str(elems[0])


def make_const(val: Any) -> Callable:
    """Return a function that simply returns the value given as input here."""

    def const(_elems):
        return val

    return const
