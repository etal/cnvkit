"""Combiner functions for Python list-like input."""
from __future__ import print_function, absolute_import, division

import pandas as pd


def first_of(elems):
    """Return the first element of the input."""
    return elems[0]


def last_of(elems):
    """Return the last element of the input."""
    return elems[-1]


def join_strings(elems, sep=','):
    """Join a Series of strings by commas."""
    # ENH if ser elements are also comma-separated, split+uniq those too
    return sep.join(pd.unique(elems))


def merge_strands(elems):
    strands = set(elems)
    if len(strands) > 1:
        return '.'
    return elems[0]


def make_const(val):
    def const(elems):
        return val
    return const

