#!/usr/bin/env python

"""Jenks natural breaks optimization.

See:

- https://en.wikipedia.org/wiki/Jenks_natural_breaks_optimization
- https://www.macwright.org/2013/02/18/literate-jenks.html
- https://www.macwright.org/simple-statistics/docs/simple_statistics.html
"""
from __future__ import division, print_function

import numpy as np


def jenks(data, n_classes):
    """Find breakpoints using Jenks natural breaks optimization.

    NB: The `data` array is sorted in-place.
    """
    assert n_classes <= len(data)
    data.sort()
    lower_class_limits, _variance_combinations = jenks_matrices(data, n_classes)
    # print(variance_combinations)

    # Extract n_classes out of the computed matrices
    breakpoints = jenks_breaks(data, lower_class_limits, n_classes)
    return breakpoints


def jenks_matrices(data, n_classes):
    """Compute Matrices for Jenks.

    Compute the matrices required for Jenks breaks. These matrices can be used
    for any classing of data with `classes <= n_classes`.
    """
    n_rows = len(data) + 1
    n_cols = n_classes + 1
    # Optimal lower class limits
    lower_class_limits = np.zeros((n_rows, n_cols),
                                  dtype=np.float_)
    lower_class_limits[1,1:] = 1.0
    lower_class_limits[2:,1] = 1.0
    # Looks like e.g. if n_classes=3, len(data)=5:
    # 0 0 0 0
    # 0 1 1 1
    # 0 1 0 0
    # 0 1 0 0
    # 0 1 0 0
    # 0 1 0 0

    # Optimal variance combinations for all classes
    variance_combinations = np.zeros((n_rows, n_cols),
                                     dtype=np.float_)
    # variance_combinations[1,1:] = 0.0
    variance_combinations[2:,1:] = float('inf')
    # e.g.:
    # 0 0 0 0
    # 0 0 0 0
    # 0 i i i
    # 0 i i i
    # 0 i i i
    # 0 i i i

    for l in xrange(2, n_rows):
        # Estimate variance for each potential classing of the data, for
        # each potential number of classes.
        # ENH: lift these out of the loop & apply more array magic
        vals = data[np.arange(l - 1, -1, -1)]
        sums = np.cumsum(vals)
        sums_square = np.cumsum(vals ** 2)
        variances = [(sums_square[i] - (sums[i] ** 2)) / (i+1)
                     for i in range(l)]
        # Record the last variance in column 1
        variance_combinations[l][1] = variances[l-1]

        for m in xrange(1, l):
            # For each column from index 1 to just before the diagonal
            # Or, for each element up to the current l marker again
            lower_class_limit = l - m + 1  # >=2, <=l
            iv = lower_class_limit - 1  # >=1, <l
            variance = variances[m-1]
            for j in xrange(2, n_cols):
                # If adding this element to an existing class will increase
                # its variance beyond the limit, break the class at this
                # point, setting the `lower_class_limit` at this point.
                var_up_left = variance + variance_combinations[iv,j-1]
                if variance_combinations[l,j] >= var_up_left:
                    lower_class_limits[l,j] = lower_class_limit
                    variance_combinations[l,j] = var_up_left

    # Return the two matrices.
    # Only `lower_class_limits` is needed to calculate breaks, but
    # variances can be useful to evaluate goodness of fit.
    return lower_class_limits, variance_combinations


def jenks_breaks(data, lower_class_limits, n_classes):
    """Pull Breaks Values for Jenks.

    The second part of the jenks recipe: Take the calculated matrices and derive
    an array of n breaks.

    Backtracking, in DP lingo.
    """
    breakpoints = np.zeros(n_classes + 1)
    breakpoints[n_classes] = data[-1]  # Upper bound
    # Use the lower_class_limits matrix as indices into itself, iteratively
    lower_limit_idx = len(data)
    for j in xrange(n_classes, 0, -1):
        lower_limit_idx = lower_class_limits[lower_limit_idx, j] - 1
        breakpoints[j - 1] = data[lower_limit_idx]
    return breakpoints


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        a = np.loadtxt(sys.argv[1])
    else:
        np.random.seed(0x5EED)
        a  = np.concatenate([
            np.random.randn(600)*1.5 - 4,
            np.random.randn(400)*.5 - 1,
                             np.random.randn(900) + 2])
    print("Data range:", a.min(), a.max())
    bp = jenks(a, 3)
    print("Breakpoints:", *bp)
    lower = a[a < bp[1]]
    middle = a[(bp[1] <= a) & (a <= bp[2])]
    upper = a[bp[2] < a]
    print("Set sizes:", len(lower), len(middle), len(upper))
    print("Total:", sum([len(lower), len(middle), len(upper)]))

