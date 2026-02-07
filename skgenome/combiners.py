"""Combiner functions for Python list-like input."""

from typing import TYPE_CHECKING, Any, Optional
from collections.abc import Callable
from collections.abc import Iterable, Sequence

import pandas as pd

if TYPE_CHECKING:
    from pandas.core.frame import DataFrame


def get_combiners(
    table: pd.DataFrame,
    stranded: bool = False,
    combine: Optional[dict[str, Callable]] = None,
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


def join_strings(elems: Iterable, sep: str = ",") -> str:
    """Join a Series of strings by commas."""
    # ENH if elements are also comma-separated, split+uniq those too
    return sep.join(pd.unique(elems))


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
