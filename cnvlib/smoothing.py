"""Signal smoothing functions."""
from __future__ import absolute_import, division

import math

import numpy as np
import pandas as pd

from . import metrics


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
    # Pad the edges of the original array with mirror copies
    signal = pd.Series(np.concatenate((x[wing-1::-1], x, x[:-wing-1:-1])))
    return x, wing, signal


def rolling_median(x, width):
    """Rolling median with mirrored edges."""
    x, wing, signal = check_inputs(x, width)
    rolled = signal.rolling(2 * wing + 1, 1, center=True).median()
    # if rolled.hasnans:
    #     rolled = rolled.interpolate()
    return np.asfarray(rolled[wing:-wing])


def rolling_quantile(x, width, quantile):
    """Rolling quantile (0--1) with mirrored edges."""
    x, wing, signal = check_inputs(x, width)
    rolled = signal.rolling(2 * wing + 1, 2, center=True).quantile(quantile)
    return np.asfarray(rolled[wing:-wing])


def rolling_std(x, width):
    """Rolling quantile (0--1) with mirrored edges."""
    x, wing, signal = check_inputs(x, width)
    rolled = signal.rolling(2 * wing + 1, 2, center=True).std()
    return np.asfarray(rolled[wing:-wing])


def smoothed(x, width, do_fit_edges=False):
    """Smooth the values in `x` with the Kaiser windowed filter.

    See: https://en.wikipedia.org/wiki/Kaiser_window

    Parameters
    ----------
    x : array-like
        1-dimensional numeric data set.
    width : float
        Fraction of x's total length to include in the rolling window (i.e. the
        proportional window width), or the integer size of the window.
    """
    x, wing, signal = check_inputs(x, width)
    # Apply signal smoothing
    window = np.kaiser(2 * wing + 1, 14)
    y = np.convolve(window / window.sum(), signal, mode='same')
    # Chop off the ends of the result so it has the original size
    y = y[wing:-wing]
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


# Outlier detection

def outlier_iqr(a, c=3.0):
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

    X_i is an outlier if::

         | X_i - M |
        _____________  > K ~= 2.24

         MAD / 0.6745

    where $K = sqrt( X^2_{0.975,1} )$,
    the square root of the 0.975 quantile of a chi-squared distribution with 1
    degree of freedom.

    This is a very robust rule with the highest possible breakdown point of 0.5.

    Returns
    -------
    np.array
        A boolean array of the same size as `a`, where outlier indices are True.

    References
    ----------
    - Davies & Gather (1993) The Identification of Multiple Outliers.
    - Rand R. Wilcox (2012) Introduction to robust estimation and hypothesis
      testing. Ch.3: Estimating measures of location and scale.
    """
    K = 2.24

    a = np.asarray(a)
    dists = np.abs(a - np.median(a))
    mad = metrics.median_absolute_deviation(a)
    return (dists / mad) > K


def rolling_outlier_iqr(x, width, c=3.0):
    """Detect outliers as a multiple of the IQR from the median.

    By convention, "outliers" are points more than 1.5 * IQR from the median (~2
    SD if values are normally distributed), and "extremes" or extreme outliers
    are those more than 3.0 * IQR (~4 SD).
    """
    if len(x) <= width:
        return np.zeros(len(x), dtype=np.bool_)
    dists = x - smoothed(x, width)
    q_hi = rolling_quantile(dists, width, .75)
    q_lo = rolling_quantile(dists, width, .25)
    iqr = q_hi - q_lo
    outliers = (np.abs(dists) > iqr * c)
    return outliers


def rolling_outlier_quantile(x, width, q, m):
    """Detect outliers by multiples of a quantile in a window.

    Outliers are the array elements outside `m` times the `q`'th
    quantile of deviations from the smoothed trend line, as calculated from
    the trend line residuals. (For example, take the magnitude of the 95th
    quantile times 5, and mark any elements greater than that value as
    outliers.)

    This is the smoothing method used in BIC-seq (doi:10.1073/pnas.1110574108)
    with the parameters width=200, q=.95, m=5 for WGS.

    Returns
    -------
    np.array
        A boolean array of the same size as `x`, where outlier indices are True.
    """
    if len(x) <= width:
        return np.zeros(len(x), dtype=np.bool_)
    dists = np.abs(x - smoothed(x, width))
    quants = rolling_quantile(dists, width, q)
    outliers = (dists > quants * m)
    return outliers


def rolling_outlier_std(x, width, stdevs):
    """Detect outliers by stdev within a rolling window.

    Outliers are the array elements outside `stdevs` standard deviations from
    the smoothed trend line, as calculated from the trend line residuals.

    Returns
    -------
    np.array
        A boolean array of the same size as `x`, where outlier indices are True.
    """
    if len(x) <= width:
        return np.zeros(len(x), dtype=np.bool_)
    dists = x - smoothed(x, width)
    x_std = rolling_std(dists, width)
    outliers = (np.abs(dists) > x_std * stdevs)
    return outliers
