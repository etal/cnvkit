"""Signal-smoothing functions."""
import logging
import math

import numpy as np
import pandas as pd
from scipy.signal import savgol_coeffs, savgol_filter

from . import descriptives


def check_inputs(x, width, as_series=True, weights=None):
    """Transform width into a half-window size.

    `width` is either a fraction of the length of `x` or an integer size of the
    whole window. The output half-window size is truncated to the length of `x`
    if needed.
    """
    x = np.asfarray(x)
    wing = _width2wing(width, x)
    signal = _pad_array(x, wing)
    if as_series:
        signal = pd.Series(signal)
    if weights is None:
        return x, wing, signal

    weights = _pad_array(weights, wing)
    # Linearly roll-off weights in mirrored wings
    weights[:wing] *= np.linspace(1/wing, 1, wing)
    weights[-wing:] *= np.linspace(1, 1/wing, wing)
    if as_series:
        weights = pd.Series(weights)
    return x, wing, signal, weights


def _width2wing(width, x, min_wing=3):
    """Convert a fractional or absolute width to integer half-width ("wing").
    """
    if 0 < width < 1:
        wing = int(math.ceil(len(x) * width * 0.5))
    elif width >= 2 and int(width) == width:
        # Ensure window width <= len(x) to avoid TypeError
        width = min(width, len(x) - 1)
        wing = int(width // 2)
    else:
        raise ValueError("width must be either a fraction between 0 and 1 "
                         "or an integer greater than 1 (got %s)" % width)
    wing = max(wing, min_wing)
    wing = min(wing, len(x) - 1)
    assert wing >= 1, "Wing must be at least 1 (got %s)" % wing
    return wing


def _pad_array(x, wing):
    """Pad the edges of the input array with mirror copies."""
    return np.concatenate((x[wing-1::-1],
                           x,
                           x[:-wing-1:-1]))


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


def convolve_weighted(window, signal, weights, n_iter=1):
    """Convolve a weighted window over a weighted signal array.

    Source: https://stackoverflow.com/a/46232913/10049
    """
    assert len(weights) == len(signal), (
        "len(weights) = %d, len(signal) = %d, window_size = %s"
        % (len(weights), len(signal), len(window)))
    y, w = signal, weights
    window /= window.sum()
    for _i in range(n_iter):
        logging.debug("Iteration %d: len(y)=%d, len(w)=%d",
        _i, len(y), len(w))
        D = np.convolve(w * y, window, mode='same')
        N = np.convolve(w, window, mode='same')
        y = D / N
        # Update weights to account for the smoothing
        w = np.convolve(window, w, mode='same')
    return y, w


def convolve_unweighted(window, signal, wing, n_iter=1):
    """Convolve a weighted window over array `signal`.

    Input array is assumed padded by `_pad_array`; output has padding removed.
    """
    window /= window.sum()
    y = signal
    for _i in range(n_iter):
        y = np.convolve(window, y, mode='same')
    # Chop off the ends of the result so it has the original size
    y = y[wing:-wing]
    return y


def guess_window_size(x, weights=None):
    """Choose a reasonable window size given the signal.

    Inspired by Silverman's rule: bandwidth is proportional to signal's standard
    deviation and the length of the signal ^ 4/5.
    """
    if weights is None:
        sd = descriptives.biweight_midvariance(x)
    else:
        sd = descriptives.weighted_std(x, weights)
    width = 4 * sd * len(x) ** (4/5)
    width = max(3, int(round(width)))
    width = min(len(x), width)
    return width


def kaiser(x, width=None, weights=None, do_fit_edges=False):
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
    if len(x) < 2:
        return x
    if width is None:
        width = guess_window_size(x, weights)
    x, wing, *padded = check_inputs(x, width, False, weights)
    # Apply signal smoothing
    window = np.kaiser(2 * wing + 1, 14)
    if weights is None:
        signal, = padded
        y = convolve_unweighted(window, signal, wing)
    else:
        signal, weights = padded
        y, _w = convolve_weighted(window, signal, weights)
    if do_fit_edges:
        _fit_edges(x, y, wing)  # In-place
    return y


def savgol(x, total_width=None, weights=None,
           window_width=7, order=3, n_iter=1):
    """Savitzky-Golay smoothing.

    Fitted polynomial order is typically much less than half the window width.

    `total_width` overrides `n_iter`.
    """
    if len(x) < 2:
        return x

    # If the effective (total) window width is not specified explicitly, compute it.
    if total_width is None:
        total_width = n_iter * window_width

    # Pad the signal.
    if weights is None:
        x, total_wing, signal = check_inputs(x, total_width, False)
    else:
        x, total_wing, signal, weights = check_inputs(x, total_width, False, weights)

    # If the signal is short, the effective window length originally requested may not be possible. Because of this, we
    # recalculate it given the actual wing length obtained.
    total_width = 2 * total_wing + 1

    # In case the signal is *very* short, the smoothing parameters will have to be adjusted as well.
    window_width = min(window_width, total_width)
    order = min(order, window_width // 2)

    # Given the adjusted window widths (one-iteration and total), calculate the number of iterations we can do.
    n_iter = max(1, min(1000, total_width // window_width))

    # Apply signal smoothing.
    logging.debug('Smoothing in {} iterations with window size {} and order {} for effective bandwidth {}',
                  n_iter, window_width, order, total_width)
    if weights is None:
        y = signal
        for _i in range(n_iter):
            y = savgol_filter(y, window_width, order, mode='interp')
        # y = convolve_unweighted(window, signal, wing)
    else:
        # TODO fit edges here, too
        window = savgol_coeffs(window_width, order)
        y, w = convolve_weighted(window, signal, weights, n_iter)
    # Safety
    bad_idx = (y > x.max()) | (y < x.min())
    if bad_idx.any():
        logging.warning("Smoothing overshot at {} / {} indices: "
                        "({}, {}) vs. original ({}, {})"
                        .format(bad_idx.sum(), len(bad_idx),
                                y.min(), y.max(),
                                x.min(), x.max()))
    return y[total_wing:-total_wing]


def _fit_edges(x, y, wing, polyorder=3):
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
    iqr = descriptives.interquartile_range(a)
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
    mad = descriptives.median_absolute_deviation(a)
    return (dists / mad) > K


def rolling_outlier_iqr(x, width, c=3.0):
    """Detect outliers as a multiple of the IQR from the median.

    By convention, "outliers" are points more than 1.5 * IQR from the median (~2
    SD if values are normally distributed), and "extremes" or extreme outliers
    are those more than 3.0 * IQR (~4 SD).
    """
    if len(x) <= width:
        return np.zeros(len(x), dtype=np.bool_)
    dists = x - savgol(x, width)
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
    dists = np.abs(x - savgol(x, width))
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
    dists = x - savgol(x, width)
    x_std = rolling_std(dists, width)
    outliers = (np.abs(dists) > x_std * stdevs)
    return outliers
