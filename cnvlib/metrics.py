"""Robust estimators of central tendency and scale.

For use in evaluating performance of copy number estimation.

See:
    https://en.wikipedia.org/wiki/Robust_measures_of_scale
    https://astropy.readthedocs.io/en/latest/_modules/astropy/stats/funcs.html

"""
from __future__ import division

import logging

import numpy as np
from scipy import stats


def ests_of_scale(deviations):
    """Estimators of scale: standard deviation, MAD, biweight midvariance.

    Calculates all of these values for an array of deviations and returns them
    as a tuple.
    """
    std = np.std(deviations, dtype=np.float64)
    mad = median_absolute_deviation(deviations)
    iqr = interquartile_range(deviations)
    biw = biweight_midvariance(deviations)
    return (std, mad, iqr, biw)


# M-estimators of central location

def biweight_location(a, initial=None, c=6.0, epsilon=1e-4):
    """Compute the biweight location for an array.

    The biweight is a robust statistic for determining the central location of a
    distribution.
    """
    a = np.asarray(a)
    if initial is None:
        initial = np.median(a)
    # Weight the observations by distance from initial estimate
    d = a - initial
    w = d / max(c * median_absolute_deviation(a), epsilon)
    w = (1 - w**2)**2
    # Omit the outlier points
    mask = (w < 1)
    weightsum = w[mask].sum()
    if weightsum == 0:
        # Insufficient variation to improve the initial estimate
        return initial
    return initial + (d[mask] * w[mask]).sum() / weightsum


def modal_location(arr):
    """Return the modal value of an array's values.

    The "mode" is the location of peak density among the values, estimated using
    a Gaussian kernel density estimator.

    `arr` is a 1-D array of floating-point values, e.g. bin log2 ratio values.
    """
    sarr = np.sort(arr)
    kde = stats.gaussian_kde(sarr)
    y = kde.evaluate(sarr)
    peak = sarr[y.argmax()]
    return peak


def segment_mean(cnarr, skip_low=False):
    """Weighted average of bin log2 values."""
    if skip_low:
        cnarr = cnarr.drop_low_coverage()
    if len(cnarr) == 0:
        return np.nan
    if 'weight' in cnarr:
        return np.average(cnarr['log2'], weights=cnarr['weight'])
    return cnarr['log2'].mean()


# Estimators of scale

def biweight_midvariance(a, initial=None, c=9.0, epsilon=1e-4):
    """Compute the biweight midvariance for an array.

    The biweight midvariance is a robust statistic for determining the
    midvariance (i.e. the standard deviation) of a distribution.

    See:
    https://en.wikipedia.org/wiki/Robust_measures_of_scale#The_biweight_midvariance
    https://astropy.readthedocs.io/en/latest/_modules/astropy/stats/funcs.html
    """
    a = np.asarray(a)
    if initial is None:
        initial = np.median(a)
    # Difference of observations from initial estimate
    d = a - initial
    # Weighting (avoid dividing by zero)
    w = d / max(c * median_absolute_deviation(a), epsilon)
    w = w**2
    # Omit the outlier points
    mask = np.abs(w) < 1
    n = mask.sum()
    return (n**0.5 * (d[mask] * d[mask] * (1 - w[mask])**4).sum()**0.5
            / np.abs(((1 - w[mask]) * (1 - 5 * w[mask])).sum()))


def interquartile_range(a):
    """Compute the difference between the array's first and third quartiles."""
    a = np.asarray(a)
    return np.percentile(a, 75) - np.percentile(a, 25)


def median_absolute_deviation(a, scale_to_sd=True):
    """Compute the median absolute deviation (MAD) of array elements.

    The MAD is defined as: ``median(abs(a - median(a)))``.

    See: https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    a = np.asarray(a)
    a_median = np.median(a)
    mad = np.median(np.abs(a - a_median))
    if scale_to_sd:
        mad *= 1.4826
    return mad


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
    a = np.asarray(a)

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


# Intervals

def confidence_interval_bootstrap(bins, alpha, bootstraps=100, smoothed=True):
    """Confidence interval for segment mean log2 value, estimated by bootstrap."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1; got %s" % alpha)
    if bootstraps <= 2 / alpha:
        new_boots = int(2 / alpha + 0.5)
        logging.warn("%d bootstraps not enough to estimate CI alpha level %f; "
                     "increasing to %d", bootstraps, alpha, new_boots)
        bootstraps = new_boots
    # Bootstrap for CI
    k = len(bins)
    rand_indices = np.random.randint(0, k, (bootstraps, k))
    samples = [bins.data.take(idx) for idx in rand_indices]
    if smoothed:
        # Essentially, resample from a kernel density estimate of the data
        # instead of the original data.
        # Silverman's Rule for KDE bandwidth (roughly):
        # std = biweight_midvariance(bins['log2'])
        std = interquartile_range(bins['log2']) / 1.34
        bw = std * (k*3/4) ** (-1/5)
        samples = [samp.assign(log2=lambda x: x['log2'] + bw * np.random.randn(k))
                   for samp in samples]
    # Recalculate segment means
    bootstrap_dist = np.array([segment_mean(samp) for samp in samples])
    alphas = np.array([alpha / 2, 1 - alpha / 2])
    if not smoothed:
        # BCa correction (Efron 1987, "Better Bootstrap Confidence Intervals")
        # http://www.tandfonline.com/doi/abs/10.1080/01621459.1987.10478410
        # Ported from R package "bootstrap" function "bcanon"
        orig_mean = segment_mean(bins)
        z0 = stats.norm.ppf((bootstrap_dist < orig_mean).sum() / bootstraps)
        zalpha = stats.norm.ppf(alphas)
        # Jackknife influence values
        u = np.array([segment_mean(bins.concat([bins[:i], bins[i+1:]]))
                      for i in range(len(bins))])
        uu = u.mean() - u
        acc = (u**3).sum() / (6 * (uu**2).sum()**1.5)
        alphas = stats.norm.cdf(z0 + (z0 + zalpha)
                                        / (1 - acc * (z0 + zalpha)))
    ci = np.percentile(bootstrap_dist, 100 * alphas)
    return ci


def prediction_interval(bins, alpha):
    """Prediction interval, estimated by percentiles."""
    pct_lo = 100 * alpha / 2
    pct_hi = 100 * (1 - alpha / 2)
    # ENH: weighted percentile
    return np.percentile(bins['log2'], [pct_lo, pct_hi])

