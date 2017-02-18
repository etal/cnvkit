"""Operations on chromosome/contig/sequence names."""
from __future__ import absolute_import, division, print_function

from itertools import takewhile

import numpy as np
import pandas as pd


def detect_big_chroms(sizes):
    """Determine the number of "big" chromosomes from their lengths.

    In the human genome, this returns 24, where the canonical chromosomes 1-22,
    X, and Y are considered "big", while mitochrondria and the alternative
    contigs are not. This allows us to exclude the non-canonical chromosomes
    from an analysis where they're not relevant.

    Returns
    -------
    n_big : int
        Number of "big" chromosomes in the genome.
    thresh : int
        Length of the smallest "big" chromosomes.
    """
    sizes = pd.Series(sizes).sort_values(ascending=False)
    reldiff = sizes.diff().abs().values[1:] / sizes.values[:-1]
    changepoints = np.nonzero(reldiff > .5)[0]
    if changepoints.any():
        n_big = changepoints[0] + 1
        thresh = sizes.iat[n_big - 1]
    else:
        n_big = len(sizes)
        thresh = sizes[-1]
    return n_big, thresh


def sorter_chrom(label):
    """Create a sorting key from chromosome label.

    Sort by integers first, then letters or strings. The prefix "chr"
    (case-insensitive), if present, is stripped automatically for sorting.

    E.g. chr1 < chr2 < chr10 < chrX < chrY < chrM
    """
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key
