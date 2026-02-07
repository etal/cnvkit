"""DataFrame-level subtraction operations.

Subtract one set of regions from another, returning the one-way difference.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""

from __future__ import annotations

from typing import TYPE_CHECKING

import bioframe

if TYPE_CHECKING:
    import pandas as pd

_BF_COLS = ("chromosome", "start", "end")


def subtract(table: pd.DataFrame, other: pd.DataFrame) -> pd.DataFrame:
    """Subtract one set of regions from another, returning the one-way difference."""
    if not len(other):
        return table
    return bioframe.subtract(table, other, cols1=_BF_COLS, cols2=_BF_COLS)
