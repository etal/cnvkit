"""Robust metrics to evaluate performance of copy number estimates.
"""
from __future__ import absolute_import, division, print_function
from builtins import map, range, zip

import logging

import numpy as np
# import pandas as pd
from scipy import stats

from . import descriptives


def do_segmetrics(cnarr, segarr, location_stats=(), spread_stats=(),
                  interval_stats=(), alpha=.05, bootstraps=100):
    """Compute segment-level metrics from bin-level log2 ratios."""
    # Silence sem's "Degrees of freedom <= 0 for slice"; NaN is OK
    import warnings
    warnings.simplefilter('ignore', RuntimeWarning)

    stat_funcs = {
        'mean': np.mean,
        'median': np.median,
        'mode': descriptives.modal_location,

        'stdev': np.std,
        'mad':  descriptives.median_absolute_deviation,
        'mse':  descriptives.mean_squared_error,
        'iqr':  descriptives.interquartile_range,
        'bivar': descriptives.biweight_midvariance,
        'sem': stats.sem,

        'ci': make_ci_func(alpha, bootstraps),
        'pi': make_pi_func(alpha),
    }

    bins_log2s = list(cnarr.iter_ranges_of(segarr, 'log2', 'outer', True))
    segarr = segarr.copy()
    if location_stats:
        # Measures of location
        for statname in location_stats:
            func = stat_funcs[statname]
            segarr[statname] = np.fromiter(map(func, bins_log2s),
                                           np.float_, len(segarr))
    # Measures of spread
    if spread_stats:
        deviations = (bl - sl for bl, sl in zip(bins_log2s, segarr['log2']))
        if len(spread_stats) > 1:
            deviations = list(deviations)
        for statname in spread_stats:
            func = stat_funcs[statname]
            segarr[statname] = np.fromiter(map(func, deviations),
                                           np.float_, len(segarr))
    # Interval calculations
    weights = cnarr['weight']
    if 'ci' in interval_stats:
        segarr['ci_lo'], segarr['ci_hi'] = calc_intervals(bins_log2s, weights,
                                                          stat_funcs['ci'])
    if 'pi' in interval_stats:
        segarr['pi_lo'], segarr['pi_hi'] = calc_intervals(bins_log2s, weights,
                                                          stat_funcs['pi'])

    return segarr


def make_ci_func(alpha, bootstraps):
    def ci_func(ser, wt):
        return confidence_interval_bootstrap(ser, wt, alpha, bootstraps)
    return ci_func


def make_pi_func(alpha):
    """Prediction interval, estimated by percentiles."""
    # ENH: weighted percentile
    pct_lo = 100 * alpha / 2
    pct_hi = 100 * (1 - alpha / 2)
    def pi_func(ser, _w):
        return np.percentile(ser, [pct_lo, pct_hi])
    return pi_func


def calc_intervals(bins_log2s, weights, func):
    """Compute a stat that yields intervals (low & high values)"""
    out_vals_lo =  np.repeat(np.nan, len(bins_log2s))
    out_vals_hi = np.repeat(np.nan, len(bins_log2s))
    for i, ser in enumerate(bins_log2s):
        if len(ser):
            wt = weights[ser.index]
            assert (wt.index == ser.index).all()
            out_vals_lo[i], out_vals_hi[i] = func(ser.values, wt.values)
    return out_vals_lo, out_vals_hi


def confidence_interval_bootstrap(values, weights, alpha, bootstraps=100, smoothed=True):
    """Confidence interval for segment mean log2 value, estimated by bootstrap."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1; got %s" % alpha)
    if bootstraps <= 2 / alpha:
        new_boots = int(np.ceil(2 / alpha))
        logging.warning("%d bootstraps not enough to estimate CI alpha level "
                        "%f; increasing to %d", bootstraps, alpha, new_boots)
        bootstraps = new_boots
    # Bootstrap for CI
    k = len(values)
    if k < 2:
        return np.repeat(values[0], 2)

    np.random.seed(0xA5EED)
    rand_indices = np.random.randint(0, k, size=(bootstraps, k))
    samples = ((np.take(values, idx), np.take(weights, idx))
               for idx in rand_indices)
    if smoothed:
        # samples = _smooth_samples(values, samples, alpha)
        pass
    # Recalculate segment means
    seg_means = (np.average(val, weights=wt)
                 for val, wt in samples)
    bootstrap_dist = np.fromiter(seg_means, np.float_, bootstraps)
    alphas = np.array([alpha / 2, 1 - alpha / 2])
    if not smoothed:
        # alphas = _bca_correct_alpha(values, weights, bootstrap_dist, alphas)
        pass
    ci = np.percentile(bootstrap_dist, list(100 * alphas))
    return ci


def _smooth_samples(values, samples, alpha):
    k = len(values)
    # Essentially, resample from a kernel density estimate of the data
    # instead of the original data.
    # Estimate KDE bandwidth (Polansky 1995)
    resids = values - values.mean()
    s_hat = 1/k * (resids**2).sum()  # sigma^2 = E[X-theta]^2
    y_hat = 1/k * abs((resids**3).sum())  # gamma = E[X-theta]^3
    z = stats.norm.ppf(alpha / 2)  # or alpha?
    bw = k**(-1/4) * np.sqrt(y_hat*(z**2 + 2) / (3*s_hat*z))
    # NB: or, Silverman's Rule for KDE bandwidth (roughly):
    # std = interquartile_range(values) / 1.34
    # bw = std * (k*3/4) ** (-1/5)
    if bw > 0:
        # Unpack the (log-ratio, weight) tuple to retain weights
        samples = [(v + bw * np.random.randn(k), w)
                   for v, w in samples]
        logging.debug("Smoothing worked for this segment (bw=%s)", bw)
    else:
        logging.debug("Smoothing not needed for this segment (bw=%s)", bw)
    return samples


def _bca_correct_alpha(values, weights, bootstrap_dist, alphas):
    # BCa correction (Efron 1987, "Better Bootstrap Confidence Intervals")
    # http://www.tandfonline.com/doi/abs/10.1080/01621459.1987.10478410
    # Ported from R package "bootstrap" function "bcanon"
    n_boots = len(bootstrap_dist)
    orig_mean = np.average(values, weights=weights)
    logging.warning("boot samples less: %s / %s",
                    (bootstrap_dist < orig_mean).sum(),
                    n_boots)
    n_boots_below = (bootstrap_dist < orig_mean).sum()
    if n_boots_below == 0:
        logging.warning("boots mean %s, orig mean %s",
                        bootstrap_dist.mean(), orig_mean)
    else:
        logging.warning("boot samples less: %s / %s",
                        n_boots_below, n_boots)
    z0 = stats.norm.ppf((bootstrap_dist < orig_mean).sum() / n_boots)
    zalpha = stats.norm.ppf(alphas)
    # Jackknife influence values
    u = np.array([np.average(np.concatenate([values[:i], values[i+1:]]),
                             weights=np.concatenate([weights[:i],
                                                     weights[i+1:]]))
                  for i in range(len(values))])
    uu = u.mean() - u
    acc = (u**3).sum() / (6 * (uu**2).sum()**1.5)
    alphas = stats.norm.cdf(z0 + (z0 + zalpha)
                                    / (1 - acc * (z0 + zalpha)))
    logging.warning("New alphas: %s -- via z0=%s, za=%s, acc=%s",
                    alphas, z0, zalpha, acc)
    if not 0 < alphas[0] < 1 and 0 < alphas[1] < 1:
        raise ValueError("CI alphas should be in (0,1); got %s" % alphas)
    return alphas


def segment_mean(cnarr, skip_low=False):
    """Weighted average of bin log2 values."""
    if skip_low:
        cnarr = cnarr.drop_low_coverage()
    if len(cnarr) == 0:
        return np.nan
    if 'weight' in cnarr:
        return np.average(cnarr['log2'], weights=cnarr['weight'])
    return cnarr['log2'].mean()
