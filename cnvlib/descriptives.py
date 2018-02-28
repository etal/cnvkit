"""Robust estimators of central tendency and scale.

See:
    https://en.wikipedia.org/wiki/Robust_measures_of_scale
    https://astropy.readthedocs.io/en/latest/_modules/astropy/stats/funcs.html

"""
from __future__ import division
import sys
from functools import wraps

import numpy as np
from scipy import stats


# Decorators to coerce input and short-circuit trivial cases

def on_array(default=None):
    """Ensure `a` is a numpy array with no missing/NaN values."""
    def outer(f):
        @wraps(f)
        def wrapper(a, **kwargs):
            a = np.asfarray(a)
            a = a[~np.isnan(a)]
            if not len(a):
                return np.nan
            if len(a) == 1:
                if default is None:
                    return a[0]
                return default
            return f(a, **kwargs)
        return wrapper
    return outer


def on_weighted_array(default=None):
    """Ensure `a` and `w` are equal-length numpy arrays with no NaN values.

    For weighted descriptives -- `a` is the array of values, `w` is weights.

    1. Drop any cells in `a` that are NaN from both `a` and `w`
    2. Replace any remaining NaN cells in `w` with 0.
    """
    def outer(f):
        @wraps(f)
        def wrapper(a, w, **kwargs):
            if len(a) != len(w):
                raise ValueError("Unequal array lengths: a=%d, w=%d"
                                % (len(a), len(w)))
            if not len(a):
                return np.nan
            a = np.asfarray(a)
            w = np.asfarray(w)
            # Drop a's NaN indices from both arrays
            a_nan = np.isnan(a)
            if a_nan.any():
                a = a[~a_nan]
                if not len(a):
                    return np.nan
                w = w[~a_nan]
            if len(a) == 1:
                if default is None:
                    return a[0]
                return default
            # Fill w's NaN indices
            w_nan = np.isnan(w)
            if w_nan.any():
                w[w_nan] = 0.0
            return f(a, w, **kwargs)
        return wrapper
    return outer


# M-estimators of central location

@on_array()
def biweight_location(a, initial=None, c=6.0, epsilon=1e-3, max_iter=5):
    """Compute the biweight location for an array.

    The biweight is a robust statistic for estimating the central location of a
    distribution.
    """
    def biloc_iter(a, initial):
        # Weight the observations by distance from initial estimate
        d = a - initial
        mad = np.median(np.abs(d))
        w = d / max(c * mad, epsilon)
        w = (1 - w**2)**2
        # Omit the outlier points
        mask = (w < 1)
        weightsum = w[mask].sum()
        if weightsum == 0:
            # Insufficient variation to improve the initial estimate
            return initial
        return initial + (d[mask] * w[mask]).sum() / weightsum

    if initial is None:
        initial = np.median(a)
    for _i in range(max_iter):
        result = biloc_iter(a, initial)
        if abs(result - initial) <= epsilon:
            break
        initial = result
    return result


@on_array()
def modal_location(a):
    """Return the modal value of an array's values.

    The "mode" is the location of peak density among the values, estimated using
    a Gaussian kernel density estimator.

    Parameters
    ----------
    a : np.array
        A 1-D array of floating-point values, e.g. bin log2 ratio values.
    """
    sarr = np.sort(a)
    kde = stats.gaussian_kde(sarr)
    y = kde.evaluate(sarr)
    peak = sarr[y.argmax()]
    return peak


@on_weighted_array()
def weighted_median(a, weights):
    """Weighted median of a 1-D numeric array."""
    order = a.argsort()
    a = a[order]
    weights = weights[order]
    midpoint = 0.5 * weights.sum()
    if (weights > midpoint).any():
        # Any point with the majority of total weight must be the median
        return a[weights.argmax()]
    cumulative_weight = weights.cumsum()
    midpoint_idx = cumulative_weight.searchsorted(midpoint)
    if (midpoint_idx > 0 and
        cumulative_weight[midpoint_idx-1] - midpoint < sys.float_info.epsilon):
        # Midpoint of 2 array values
        return a[midpoint_idx-1 : midpoint_idx+1].mean()
    return a[midpoint_idx]


# Estimators of scale

@on_array(0)
def biweight_midvariance(a, initial=None, c=9.0, epsilon=1e-3):
    """Compute the biweight midvariance for an array.

    The biweight midvariance is a robust statistic for determining the
    midvariance (i.e. the standard deviation) of a distribution.

    See:

    - https://en.wikipedia.org/wiki/Robust_measures_of_scale#The_biweight_midvariance
    - https://astropy.readthedocs.io/en/latest/_modules/astropy/stats/funcs.html
    """
    if initial is None:
        initial = biweight_location(a)
    # Difference of observations from initial location estimate
    d = a - initial
    # Weighting (avoid dividing by zero)
    mad = np.median(np.abs(d))
    w = d / max(c * mad, epsilon)
    # Omit the outlier points
    mask = np.abs(w) < 1
    if w[mask].sum() == 0:
        # Insufficient variation to improve on MAD
        return mad * 1.4826
    n = mask.sum()
    d_ = d[mask]
    w_ = (w**2)[mask]
    return np.sqrt((n * (d_**2 * (1 - w_)**4).sum())
                   / (((1 - w_) * (1 - 5 * w_)).sum()**2))


@on_array(0)
def gapper_scale(a):
    """Scale estimator based on gaps between order statistics.

    See:

    - Wainer & Thissen (1976)
    - Beers, Flynn, and Gebhardt (1990)
    """
    gaps = np.diff(np.sort(a))
    n = len(a)
    idx = np.arange(1, n)
    weights = idx * (n - idx)
    return (gaps * weights).sum() * np.sqrt(np.pi) / (n * (n - 1))


@on_array(0)
def interquartile_range(a):
    """Compute the difference between the array's first and third quartiles."""
    return np.percentile(a, 75) - np.percentile(a, 25)


@on_array(0)
def median_absolute_deviation(a, scale_to_sd=True):
    """Compute the median absolute deviation (MAD) of array elements.

    The MAD is defined as: ``median(abs(a - median(a)))``.

    See: https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    a_median = np.median(a)
    mad = np.median(np.abs(a - a_median))
    if scale_to_sd:
        mad *= 1.4826
    return mad


@on_weighted_array()
def weighted_mad(a, weights, scale_to_sd=True):
    """Median absolute deviation (MAD) with weights."""
    a_median = weighted_median(a, weights)
    mad = weighted_median(np.abs(a - a_median), weights)
    if scale_to_sd:
        mad *= 1.4826
    return mad


@on_weighted_array()
def weighted_std(a, weights):
    """Standard deviation with weights."""
    mean = np.average(a, weights=weights)
    var = np.average((a - mean) ** 2, weights=weights)
    return np.sqrt(var)


@on_array(0)
def mean_squared_error(a, initial=None):
    """Mean squared error (MSE).

    By default, assume the input array `a` is the residuals/deviations/error,
    so MSE is calculated from zero. Another reference point for calculating the
    error can be specified with `initial`.
    """
    if initial is None:
        initial = a.mean()
    if initial:
        a = a - initial
    return (a ** 2).mean()


@on_array(0)
def q_n(a):
    """Rousseeuw & Croux's (1993) Q_n, an alternative to MAD.

    ``Qn := Cn first quartile of (|x_i - x_j|: i < j)``

    where Cn is a constant depending on n.

    Finite-sample correction factors must be used to calibrate the
    scale of Qn for small-to-medium-sized samples.

        n   E[Qn]
        --  -----
        10  1.392
        20  1.193
        40  1.093
        60  1.064
        80  1.048
        100 1.038
        200 1.019

    """
    # First quartile of: (|x_i - x_j|: i < j)
    vals = []
    for i, x_i in enumerate(a):
        for x_j in a[i+1:]:
            vals.append(abs(x_i - x_j))
    quartile = np.percentile(vals, 25)

    # Cn: a scaling factor determined by sample size
    n = len(a)
    if n <= 10:
        # ENH: warn when extrapolating beyond the data
        # ENH: simulate for values up to 10
        #   (unless the equation below is reliable)
        scale = 1.392
    elif 10 < n < 400:
        # I fitted the simulated values (above) to a power function in Excel:
        #   f(x) = 1.0 + 3.9559 * x ^ -1.0086
        # This should be OK for interpolation. (Does it apply generally?)
        scale = 1.0 + (4 / n)
    else:
        scale = 1.0

    return quartile / scale
