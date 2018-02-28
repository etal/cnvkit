"""Combiner functions for Python list-like input."""
from __future__ import print_function, absolute_import, division

import pandas as pd


def get_combiners(table, stranded=False, combine=None):
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
        'chromosome': first_of,
        'start': first_of,
        'end': max,
        'gene': join_strings,
        'accession': join_strings,
        'weight': sum,
        'probes': sum,
    }
    if combine:
        cmb.update(combine)
    if 'strand' not in cmb:
        cmb['strand'] = first_of if stranded else merge_strands
    return {k: v for k, v in cmb.items() if k in table.columns}


def first_of(elems):
    """Return the first element of the input."""
    return elems[0]


def last_of(elems):
    """Return the last element of the input."""
    return elems[-1]


max_of = max


def join_strings(elems, sep=','):
    """Join a Series of strings by commas."""
    # ENH if elements are also comma-separated, split+uniq those too
    return sep.join(pd.unique(elems))


def merge_strands(elems):
    strands = set(elems)
    if len(strands) > 1:
        return '.'
    return elems[0]


def make_const(val):
    def const(_elems):
        return val
    return const
