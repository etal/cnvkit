"""DataFrame-level subtraction operations.

Subtract one set of regions from another, returning the one-way difference.

The functions here operate on pandas DataFrame and Series instances, not
GenomicArray types.

"""
from __future__ import print_function, absolute_import, division
from builtins import next

import logging

import pandas as pd


def _subtract(table, other):
    if not len(other):
        return table

    other_chroms = {c: o for c, o in other.groupby(['chromosome'], sort=False)}
    done_chroms = []
    for chrom, ctable in table.groupby('chromosome', sort=False):
        if chrom in other_chroms:
            logging.info("%s: Subtracting excluded regions", chrom)
            newrows = subtract_chrom(ctable, other_chroms[chrom])
            done_chroms.append(
                pd.DataFrame.from_records(newrows, columns=table.columns))
        else:
            logging.info("%s: No excluded regions", chrom)
            done_chroms.append(ctable)

    return pd.concat(done_chroms)


def subtract_chrom(rows, other_rows):
    exclude_rows = iter(other_rows.itertuples(index=False))
    ex_start, ex_end = _next_or_inf(exclude_rows)
    for row in rows.itertuples(index=False):
        for _c, start, end in exclude_in_region(exclude_rows, row.chromosome,
                                                row.start, row.end,
                                                ex_start, ex_end):
            yield row._replace(start=start, end=end)


# TODO - rewrite non-recursively (#150)
def exclude_in_region(exclude_rows, chrom, a_start, a_end, ex_start, ex_end):
    """Take region exclusions from an iterable and apply, perhaps recursively.

    Yields
    ------
    tuple
        A pair of tuples:
            - (accessible chromosome, start, end)
            - (current exclusion start, end)
    """
    # If we've leapfrogged the excluded area, catch up
    while ex_end <= a_start:
        ex_start, ex_end = _next_or_inf(exclude_rows)
    if a_end <= ex_start:
        # Excluded area does not overlap this one
        yield (chrom, a_start, a_end)
    else:
        # Excluded area overlaps this one -> trim this region
        logging.info("\tExclusion %s:%d-%d overlaps accessible region %d-%d",
                     chrom, ex_start, ex_end, a_start, a_end)
        if ex_start <= a_start:
            if ex_end < a_end:
                # Exclusion covers this region's left (start) edge only
                for cse in exclude_in_region(exclude_rows, chrom, ex_end, a_end,
                                             ex_start, ex_end):
                    yield cse
            # Otherwise: Exclusion covers the whole region
        else:
            yield (chrom, a_start, ex_start)
            if ex_end < a_end:
                # Exclusion is in the middle of this region
                for cse in exclude_in_region(exclude_rows, chrom, ex_end,
                                             a_end, ex_start, ex_end):
                    yield cse
            # Otherwise: Exclusion covers this region's right (end) edge only


def _next_or_inf(iterable):
    try:
        ex = next(iterable)
        return ex.start, ex.end
    except StopIteration:
        return float("Inf"), float("Inf")
