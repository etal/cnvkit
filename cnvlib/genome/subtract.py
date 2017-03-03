"""DataFrame-level subtraction operations.

Subtract one set of regions from another, returning the one-way difference.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
from __future__ import print_function, absolute_import, division

import logging

import numpy as np
import pandas as pd

from .intersect import by_ranges


def subtract(table, other):
    if not len(other):
        return table
    return pd.DataFrame.from_records(_subtraction(table, other),
                                     columns=table.columns)


def _subtraction(table, other):
    for keeper, rows_to_exclude in by_ranges(other, table, 'outer', True):
        if len(rows_to_exclude):
            logging.debug(" %s:%d-%d : Subtracting %d excluded regions",
                          keeper.chromosome, keeper.start, keeper.end,
                          len(rows_to_exclude))

            keep_left = (keeper.start < rows_to_exclude.start.iat[0])
            keep_right = (keeper.end > rows_to_exclude.end.iat[-1])
            if keep_left and keep_right:
                # Keep both original edges of the source region
                # =========
                #   --  --
                starts = np.r_[keeper.start, rows_to_exclude.end.values]
                ends = np.r_[rows_to_exclude.start.values, keeper.end]
            elif keep_left:
                # Exclusion overlaps only the right side
                # =======
                #   -- ---
                starts = np.r_[keeper.start, rows_to_exclude.end.values[:-1]]
                ends = rows_to_exclude.start.values
            elif keep_right:
                # Exclusion overlaps only the left side
                #  ========
                # ---  --
                starts = rows_to_exclude.end.values
                ends = np.r_[rows_to_exclude.start.values[1:], keeper.end]
            elif len(rows_to_exclude) > 1:
                # Exclusions overlap both edges
                #  ======
                # -- -- ---
                starts = rows_to_exclude.end.values[:-1]
                ends = rows_to_exclude.start.values[1:]
            else:
                # Exclusion covers the whole region
                continue

            for start, end in zip(starts, ends):
                if end > start:
                    yield keeper._replace(start=start, end=end)
                else:
                    logging.debug("Discarding pair: (%d, %d)", start, end)

        else:
            logging.debug(" %s:%d-%d : No excluded regions",
                          keeper.chromosome, keeper.start, keeper.end)
            yield keeper
