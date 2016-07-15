"""Robust metrics to evaluate performance of copy number estimates.
"""
from __future__ import division

import logging

import numpy as np
from scipy import stats

from .descriptives import (biweight_midvariance, interquartile_range,
                           median_absolute_deviation)


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


def segment_mean(cnarr, skip_low=False):
    """Weighted average of bin log2 values."""
    if skip_low:
        cnarr = cnarr.drop_low_coverage()
    if len(cnarr) == 0:
        return np.nan
    if 'weight' in cnarr:
        return np.average(cnarr['log2'], weights=cnarr['weight'])
    return cnarr['log2'].mean()


# Intervals

def confidence_interval_bootstrap(bins, alpha, bootstraps=100, smoothed=True):
    """Confidence interval for segment mean log2 value, estimated by bootstrap."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1; got %s" % alpha)
    if bootstraps <= 2 / alpha:
        new_boots = int(np.ceil(2 / alpha))
        logging.warn("%d bootstraps not enough to estimate CI alpha level %f; "
                     "increasing to %d", bootstraps, alpha, new_boots)
        bootstraps = new_boots
    # Bootstrap for CI
    k = len(bins)
    if k < 2:
        return np.array([bins["log2"][0], bins["log2"][0]])

    rand_indices = np.random.randint(0, k, (bootstraps, k))
    samples = [bins.data.take(idx) for idx in rand_indices]
    if smoothed:
        # samples = _smooth_samples(bins, samples, alpha)
        pass
    # Recalculate segment means
    bootstrap_dist = np.array([segment_mean(samp) for samp in samples])
    alphas = np.array([alpha / 2, 1 - alpha / 2])
    if not smoothed:
        # alphas = _bca_correct_alpha(bins, bootstrap_dist, alphas)
        pass
    ci = np.percentile(bootstrap_dist, list(100 * alphas))
    return ci


def _smooth_samples(bins, samples, alpha):
    k = len(bins)
    # Essentially, resample from a kernel density estimate of the data
    # instead of the original data.
    # Estimate KDE bandwidth (Polansky 1995)
    resids = bins['log2'] - bins['log2'].mean()
    s_hat = 1/k * (resids**2).sum()  # sigma^2 = E[X-theta]^2
    y_hat = 1/k * abs((resids**3).sum())  # gamma = E[X-theta]^3
    z = stats.norm.ppf(alpha / 2)  # or alpha?
    bw = k**(-1/4) * np.sqrt(y_hat*(z**2 + 2) / (3*s_hat*z))
    # NB: or, Silverman's Rule for KDE bandwidth (roughly):
    # std = interquartile_range(bins['log2']) / 1.34
    # bw = std * (k*3/4) ** (-1/5)
    if bw > 0:
        samples = [samp.assign(log2=lambda x:
                                x['log2'] + bw * np.random.randn(k))
                    for samp in samples]
        logging.debug("Smoothing worked for this segment (bw=%s)", bw)
    else:
        logging.debug("Smoothing not needed for this segment (bw=%s)", bw)
    return samples


def _bca_correct_alpha(bins, bootstrap_dist, alphas):
    # BCa correction (Efron 1987, "Better Bootstrap Confidence Intervals")
    # http://www.tandfonline.com/doi/abs/10.1080/01621459.1987.10478410
    # Ported from R package "bootstrap" function "bcanon"
    n_boots = len(bootstrap_dist)
    orig_mean = segment_mean(bins)
    logging.warn("boot samples less: %s / %s",
                 (bootstrap_dist < orig_mean).sum(),
                 n_boots)
    n_boots_below = (bootstrap_dist < orig_mean).sum()
    if n_boots_below == 0:
        logging.warn("boots mean %s, orig mean %s",
                     bootstrap_dist.mean(), orig_mean)
    else:
        logging.warn("boot samples less: %s / %s",
                     n_boots_below, n_boots)
    z0 = stats.norm.ppf((bootstrap_dist < orig_mean).sum() / n_boots)
    zalpha = stats.norm.ppf(alphas)
    # Jackknife influence values
    u = np.array([segment_mean(bins.concat([bins[:i], bins[i+1:]]))
                    for i in range(len(bins))])
    uu = u.mean() - u
    acc = (u**3).sum() / (6 * (uu**2).sum()**1.5)
    alphas = stats.norm.cdf(z0 + (z0 + zalpha)
                                    / (1 - acc * (z0 + zalpha)))
    logging.warn("New alphas: %s -- via z0=%s, za=%s, acc=%s",
                    alphas, z0, zalpha, acc)
    if not 0 < alphas[0] < 1 and 0 < alphas[1] < 1:
        raise ValueError("CI alphas should be in (0,1); got %s" % alphas)
    return alphas


def prediction_interval(bins, alpha):
    """Prediction interval, estimated by percentiles."""
    pct_lo = 100 * alpha / 2
    pct_hi = 100 * (1 - alpha / 2)
    # ENH: weighted percentile
    return np.percentile(bins['log2'], [pct_lo, pct_hi])
