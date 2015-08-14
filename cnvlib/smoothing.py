"""Signal smoothing functions."""
from __future__ import absolute_import, division

from bisect import insort, bisect_left
from collections import deque
import math

import numpy as np

from . import core, metrics


def check_inputs(x, width):
    """Transform width into a half-window size.

    `width` is either a fraction of the length of `x` or an integer size of the
    whole window. The output half-window size is truncated to the length of `x`
    if needed.
    """
    x = np.asfarray(x)
    if 0 < width < 1:
        wing = int(math.ceil(len(x) * width * 0.5))
    elif width >= 2 and int(width) == width:
        wing = width // 2
    else:
        raise ValueError("width must be between 0 and 1 (got %s)" % width)
    if wing > len(x):
        wing = len(x) - 1
    assert wing > 0, "Wing must be greater than 0 (got %s)" % wing
    return x, wing


def rolling_median(x, width):
    """Rolling median.

    Contributed by Peter Otten to comp.lang.python.

    Source:
    https://bitbucket.org/janto/snippets/src/tip/running_median.py
    https://groups.google.com/d/msg/comp.lang.python/0OARyHF0wtA/SEs-glW4t6gJ
    """
    x, wing = check_inputs(x, width)
    # Pre-allocate the result array
    result = np.empty_like(x)
    # Keep a copy of the rolling window in original order; initially fill with
    # a mirrored copy of the first 'wing' points
    window = deque(np.concatenate((x[wing::-1], x[:wing])))
    # Also keep a sorted copy of the rolling window values
    sortwin = sorted(window)
    # Pad the right edge of the original array with a mirror copy
    signal = np.concatenate((x[wing:], x[:-wing-1:-1]))
    # Calculate the rolling median at each subsequent point
    for i, item in enumerate(signal):
        old = window.popleft()
        window.append(item)
        del sortwin[bisect_left(sortwin, old)]
        insort(sortwin, item)
        result[i] = sortwin[wing]
    return result


def smoothed(x, width, do_fit_edges=False):
    """Smooth the values in `x` with the Kaiser windowed filter.

    See: https://en.wikipedia.org/wiki/Kaiser_window

    Parameters:

    x : array-like
        1-dimensional numeric data set.
    width : float
        Fraction of x's total length to include in the rolling window (i.e. the
        proportional window width), or the integer size of the window.
    """
    x, wing = check_inputs(x, width)
    # Pad the edges with mirror-image copies of the array
    signal = np.concatenate((x[wing-1::-1], x, x[:-wing:-1]))
    # Apply signal smoothing
    window = np.kaiser(2* wing + 1, 14)
    y = np.convolve(window / window.sum(), signal, mode='same')
    # Chop off the ends of the result so it has the original size
    y = y[wing:1-wing]
    if do_fit_edges:
        fit_edges(x, y, wing)  # In-place
    return y


def fit_edges(x, y, wing, polyorder=3):
    """Apply polynomial interpolation to the edges of y, in-place.

    Calculates a polynomial fit (of order `polyorder`) of `x` within a window of
    width twice `wing`, then updates the smoothed values `y` in the half of the
    window closest to the edge.
    """
    window_length = 2 * wing + 1
    n = len(x)
    # Fit each of the two array edges (start and end)
    _fit_edge(x, y, 0, window_length, 0, wing, polyorder)
    _fit_edge(x, y, n - window_length, n, n - wing, n, polyorder)
    # TODO - fix the discontinuities at wing, n-wing


def _fit_edge(x, y, window_start, window_stop, interp_start, interp_stop,
              polyorder):
    """
    Given a 1-D array `x` and the specification of a slice of `x` from
    `window_start` to `window_stop`, create an interpolating polynomial of the
    sliced sub-array, and evaluate that polynomial from `interp_start` to
    `interp_stop`.  Put the result into the corresponding slice of `y`.
    """
    # Get the edge into a (window_length, -1) array.
    x_edge = x[window_start:window_stop]
    # Fit the edges.  poly_coeffs has shape (polyorder + 1, -1),
    # where '-1' is the same as in x_edge.
    poly_coeffs = np.polyfit(np.arange(0, window_stop - window_start),
                             x_edge, polyorder)
    # Compute the interpolated values for the edge.
    i = np.arange(interp_start - window_start, interp_stop - window_start)
    values = np.polyval(poly_coeffs, i)
    # Put the values into the appropriate slice of y.
    y[interp_start:interp_stop] = values


def smooth_genome_coverages(probes, smooth_func, width):
    """Fit a trendline through probe coverages, handling chromosome boundaries.

    Returns an array of smoothed coverage values, calculated with `smooth_func`
    and `width`, equal in length to `probes`.
    """
    # ENH: also split by centromeres (long internal gaps -- see PSCBS)
    out = {chrom: smooth_func(subprobes['log2'], width)
           for chrom, subprobes in probes.by_chromosome()}
    return np.concatenate(
        [out[chrom] for chrom in sorted(out, key=core.sorter_chrom)])


# Outlier detection

def outlier_iqr(a, c=1.5):
    """Detect outliers as a multiple of the IQR from the median.

    By convention, "outliers" are points more than 1.5 * IQR from the median,
    and "extremes" or extreme outliers are those more than 3.0 * IQR.
    """
    a = np.asarray(a)
    dists = np.abs(a - np.median(a))
    iqr = metrics.interquartile_range(a)
    return dists > (c * iqr)


def outlier_mad_median(a):
    """MAD-Median rule for detecting outliers.

    Returns: a boolean array of the same size, where outlier indices are True.

    X_i is an outlier if::

         | X_i - M |
        _____________  > K ~= 2.24

         MAD / 0.6745

    where $K = sqrt( X^2_{0.975,1} )$,
    the square root of the 0.975 quantile of a chi-squared distribution with 1
    degree of freedom.

    This is a very robust rule with the highest possible breakdown point of 0.5.

    See:

    - Davies & Gather (1993) The Identification of Multiple Outliers.
    - Rand R. Wilcox (2012) Introduction to robust estimation and hypothesis
      testing. Ch.3: Estimating measures of location and scale.
    """
    K = 2.24

    a = np.asarray(a)
    dists = np.abs(a - np.median(a))
    mad = metrics.median_absolute_deviation(a, scale_to_sd=False)
    return (dists / mad) > K
