"""DataFrame-level cutting operations.

Split genomic regions at boundaries defined by another set of regions,
similar to bedtools intersect -split.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""

from __future__ import annotations

import itertools
import numpy as np
import pandas as pd


def cut(table: pd.DataFrame, other: pd.DataFrame) -> pd.DataFrame:
    """Split intervals in ``table`` at the boundary coordinates of ``other``.

    Each interval in ``table`` is split wherever a start or end coordinate
    from ``other`` falls strictly within it. Every resulting piece inherits
    all data columns from its source row.

    Parameters
    ----------
    table : DataFrame
        Genomic intervals with at least chromosome, start, end columns.
    other : DataFrame
        Intervals whose start/end positions define the cut points.

    Returns
    -------
    DataFrame
        A copy of ``table`` with rows split at the breakpoints from ``other``.
    """
    if table.empty or other.empty:
        return table

    # Collect breakpoints per chromosome from other
    other_bp: dict[str, np.ndarray] = {}
    for chrom, grp in other.groupby("chromosome", sort=False):
        other_bp[chrom] = np.unique(
            np.concatenate([grp["start"].to_numpy(), grp["end"].to_numpy()])
        )

    result_chunks: list[pd.DataFrame] = []
    for chrom, chrom_table in table.groupby("chromosome", sort=False):
        if chrom not in other_bp:
            result_chunks.append(chrom_table)
            continue

        breakpoints = other_bp[chrom]
        cut_rows: list[tuple] = []
        for row in chrom_table.itertuples(index=False):
            # Breakpoints strictly inside (row.start, row.end)
            inner = breakpoints[(breakpoints > row.start) & (breakpoints < row.end)]
            if len(inner) == 0:
                cut_rows.append(row)
            else:
                bounds = np.concatenate([[row.start], inner, [row.end]])
                for seg_start, seg_end in itertools.pairwise(bounds):
                    cut_rows.append(
                        row._replace(start=int(seg_start), end=int(seg_end))
                    )
        result_chunks.append(
            pd.DataFrame.from_records(cut_rows, columns=chrom_table.columns)
        )

    if not result_chunks:
        return table.iloc[:0]
    return pd.concat(result_chunks, ignore_index=True)
