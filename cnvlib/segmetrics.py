"""Robust metrics to evaluate performance of copy number estimates."""

from __future__ import annotations
import logging
from typing import TYPE_CHECKING, Any

# import pandas as pd
import numpy as np
from scipy import stats

from . import descriptives

if TYPE_CHECKING:
    from collections.abc import Callable
    from collections.abc import Iterator
    from cnvlib.cnary import CopyNumArray
    from numpy import float64, ndarray
    from pandas.core.series import Series


def do_segmetrics(
    cnarr: CopyNumArray,
    segarr: CopyNumArray,
    location_stats: tuple[()] | list[str] = (),
    spread_stats: tuple[()] | list[str] = (),
    interval_stats: tuple[()] | list[str] = (),
    alpha: float = 0.05,
    bootstraps: int = 100,
    smoothed: bool | int = 10,
    skip_low: bool = False,
) -> CopyNumArray:
    """Compute segment-level metrics from bin-level log2 ratios.

    Parameters
    ----------
    cnarr : CopyNumArray
        Bin-level copy number data.
    segarr : CopyNumArray
        Segmented copy number data.
    location_stats : list of str, optional
        Location statistics to compute: 'mean', 'median', 'mode', 'p_ttest'.
        Default is empty tuple.
    spread_stats : list of str, optional
        Spread statistics to compute: 'stdev', 'mad', 'mse', 'iqr', 'bivar', 'sem'.
        Default is empty tuple.
    interval_stats : list of str, optional
        Interval statistics to compute: 'ci' (confidence interval),
        'pi' (prediction interval). Default is empty tuple.
    alpha : float, optional
        Significance level for confidence/prediction intervals. Default is 0.05.
    bootstraps : int, optional
        Number of bootstrap iterations for confidence intervals. Default is 100.
    smoothed : bool or int, optional
        Smoothed bootstrap threshold for confidence intervals. If bool: True to
        always use smoothed bootstrap, False to never use it. If int: use smoothed
        bootstrap when segment has <= this many bins. Smoothed bootstrap adds
        Gaussian noise to improve CI accuracy for small segments. BCa correction
        is applied when smoothing is not used. Default is 10.
    skip_low : bool, optional
        Skip bins with low coverage. Default is False.

    Returns
    -------
    CopyNumArray
        Segmented data with additional statistical columns.
    """
    # Silence sem's "Degrees of freedom <= 0 for slice"; NaN is OK
    import warnings

    warnings.simplefilter("ignore", RuntimeWarning)

    stat_funcs = {
        "mean": np.mean,
        "median": np.median,
        "mode": descriptives.modal_location,
        "p_ttest": lambda a: stats.ttest_1samp(a, 0.0, nan_policy="omit")[1],
        "stdev": np.std,
        "mad": descriptives.median_absolute_deviation,
        "mse": descriptives.mean_squared_error,
        "iqr": descriptives.interquartile_range,
        "bivar": descriptives.biweight_midvariance,
        "sem": stats.sem,
        "ci": make_ci_func(alpha, bootstraps, smoothed),
        "pi": make_pi_func(alpha),
    }

    if skip_low:
        cnarr = cnarr.drop_low_coverage()
    bins_log2s = list(cnarr.iter_ranges_of(segarr, "log2", "outer", True))

    segarr = segarr.copy()
    if location_stats:
        # Measures of location
        for statname in location_stats:
            func = stat_funcs[statname]
            segarr[statname] = np.fromiter(
                map(func, bins_log2s), np.float64, len(segarr)
            )
    # Measures of spread
    if spread_stats:
        deviations = (
            bl - sl for bl, sl in zip(bins_log2s, segarr["log2"], strict=True)
        )
        if len(spread_stats) > 1:
            deviations = list(deviations)  # type: ignore[assignment]
        for statname in spread_stats:
            func = stat_funcs[statname]
            segarr[statname] = np.fromiter(
                map(func, deviations), np.float64, len(segarr)
            )
    # Interval calculations
    weights = cnarr["weight"]
    if "ci" in interval_stats:
        segarr["ci_lo"], segarr["ci_hi"] = calc_intervals(
            bins_log2s, weights, stat_funcs["ci"]
        )
    if "pi" in interval_stats:
        segarr["pi_lo"], segarr["pi_hi"] = calc_intervals(
            bins_log2s, weights, stat_funcs["pi"]
        )

    return segarr


def make_ci_func(alpha: float, bootstraps: int, smoothed: bool | int) -> Callable:
    """Create a confidence interval function.

    Parameters
    ----------
    alpha : float
        Significance level for CI.
    bootstraps : int
        Number of bootstrap iterations.
    smoothed : bool or int
        If bool: True to always smooth, False to never smooth.
        If int: Threshold - smooth when n_bins <= smoothed.
    """

    def ci_func(ser, wt):
        return confidence_interval_bootstrap(ser, wt, alpha, bootstraps, smoothed)

    return ci_func


def make_pi_func(alpha: float) -> Callable:
    """Prediction interval, estimated by percentiles."""
    # ENH: weighted percentile
    pct_lo = 100 * alpha / 2
    pct_hi = 100 * (1 - alpha / 2)

    def pi_func(ser, _w):
        return np.percentile(ser, [pct_lo, pct_hi])

    return pi_func


def calc_intervals(
    bins_log2s: list[Series], weights: Series, func: Callable
) -> tuple[ndarray, ndarray]:
    """Compute a stat that yields intervals (low & high values)"""
    out_vals_lo = np.repeat(np.nan, len(bins_log2s))
    out_vals_hi = np.repeat(np.nan, len(bins_log2s))
    for i, ser in enumerate(bins_log2s):
        if len(ser):
            wt = weights[ser.index]
            assert (wt.index == ser.index).all()
            out_vals_lo[i], out_vals_hi[i] = func(ser.values, wt.values)
    return out_vals_lo, out_vals_hi


def confidence_interval_bootstrap(
    values: ndarray,
    weights: ndarray,
    alpha: float,
    bootstraps: int = 100,
    smoothed: bool | int = False,
) -> ndarray:
    """Confidence interval for segment mean log2 value, estimated by bootstrap.

    Parameters
    ----------
    values : ndarray
        Log2 ratio values.
    weights : ndarray
        Weights for each value.
    alpha : float
        Significance level for CI.
    bootstraps : int
        Number of bootstrap iterations.
    smoothed : bool or int
        If bool: True to always use smoothed bootstrap, False to never use it.
        If int: Threshold - use smoothed bootstrap when len(values) <= smoothed.
        Smoothed bootstrap adds Gaussian noise to improve CI accuracy for small
        segments. BCa correction is applied when smoothing is not used.

    Returns
    -------
    ndarray
        [ci_lo, ci_hi] confidence interval bounds.
    """
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be between 0 and 1; got {alpha}")
    if bootstraps <= 2 / alpha:
        new_boots = int(np.ceil(2 / alpha))
        logging.warning(
            "%d bootstraps not enough to estimate CI alpha level %f; increasing to %d",
            bootstraps,
            alpha,
            new_boots,
        )
        bootstraps = new_boots
    # Bootstrap for CI
    k = len(values)
    if k < 2:
        return np.repeat(values[0], 2)

    # Determine whether to use smoothed bootstrap
    if isinstance(smoothed, bool):
        use_smoothing = smoothed
    else:
        # smoothed is a threshold (integer)
        use_smoothing = k <= smoothed

    rng = np.random.default_rng(0xA5EED)
    rand_indices = rng.integers(0, k, size=(bootstraps, k))
    samples = ((np.take(values, idx), np.take(weights, idx)) for idx in rand_indices)
    if use_smoothing:
        samples = _smooth_samples_by_weight(values, samples)  # type: ignore[assignment]
    # Recalculate segment means
    seg_means = (np.average(val, weights=wt) for val, wt in samples)
    bootstrap_dist = np.fromiter(seg_means, np.float64, bootstraps)
    alphas = np.array([alpha / 2, 1 - alpha / 2])
    if not use_smoothing:
        alphas = _bca_correct_alpha(values, weights, bootstrap_dist, alphas)
    ci = np.percentile(bootstrap_dist, list(100 * alphas))
    return ci


def _smooth_samples_by_weight(
    values: ndarray, samples: Iterator[Any]
) -> list[tuple[ndarray, ndarray]]:
    """Add Gaussian noise to each bootstrap replicate.

    The result is used to compute a "smoothed bootstrap," where the added noise
    ensures that for small samples (e.g. number of bins in the segment) the
    bootstrapped CI is close to the standard error of the mean, as it should be.
    Conceptually, sample from a KDE instead of the values themselves.

    This addresses the issue that small segments (#bins < #replicates) don't
    fully represent the underlying distribution, in particular the extreme
    values, so the CI is too narrow. For single-bin segments in particular,
    the confidence interval will always have zero width unless the samples are
    smoothed.

    Standard deviation of the noise added to each bin comes from each bin's
    weight, which is an estimate of (1-variance).

    Parameters
    ----------
    values : np.ndarray
        Original log2 values within the segment.
    samples : list of np.ndarray
        Bootstrap replicates as (value_sample, weight_sample).

    Returns
    -------
    `samples` with random N(0, pop_sd) added to each value, and
    weights unchanged.
    """
    k = len(values)
    # KDE bandwidth narrows for larger sample sizes
    # Following Silverman's Rule and Polansky 1995,
    # but requiring k=1 -> bw=1 for consistency
    bw = k ** (-1 / 4)
    rng = np.random.default_rng()
    samples_list = [
        (v + (bw * np.sqrt(1 - w) * rng.standard_normal(k)), w) for v, w in samples
    ]
    return samples_list


def _bca_correct_alpha(values, weights, bootstrap_dist, alphas):
    """Bias Corrected & Accellerated (BCa) bootstrap adjustment.

    See: Efron 1987, "Better Bootstrap Confidence Intervals"
    http://www.tandfonline.com/doi/abs/10.1080/01621459.1987.10478410

    Ported from R package "bootstrap" function "bcanon".
    """
    n_boots = len(bootstrap_dist)
    orig_mean = np.average(values, weights=weights)
    n_boots_below = (bootstrap_dist < orig_mean).sum()

    # Check if proportion is too extreme for BCa
    proportion = n_boots_below / n_boots
    if proportion == 0 or proportion == 1:
        logging.warning(
            "BCa: All bootstrap samples on one side (%d/%d); using original alphas",
            n_boots_below,
            n_boots,
        )
        return alphas

    z0 = stats.norm.ppf(proportion)
    zalpha = stats.norm.ppf(alphas)

    # Jackknife influence values
    u = np.array(
        [
            np.average(
                np.concatenate([values[:i], values[i + 1 :]]),
                weights=np.concatenate([weights[:i], weights[i + 1 :]]),
            )
            for i in range(len(values))
        ]
    )
    uu = u.mean() - u
    uu_var = (uu**2).sum()

    # Check for zero variance in jackknife estimates
    if uu_var < 1e-10:
        logging.warning("BCa: Jackknife variance too small; using original alphas")
        return alphas

    acc = (uu**3).sum() / (6 * uu_var**1.5)
    denom = 1 - acc * (z0 + zalpha)

    # Check if denominator is positive
    if (denom <= 0).any():
        logging.warning(
            "BCa: Denominator non-positive (acc=%.4f); using original alphas", acc
        )
        return alphas

    new_alphas = stats.norm.cdf(z0 + (z0 + zalpha) / denom)

    # Validate new alphas
    if not (0 < new_alphas[0] < 1 and 0 < new_alphas[1] < 1):
        logging.warning(
            "BCa: Adjusted alphas %s out of range; using original alphas", new_alphas
        )
        return alphas

    logging.debug(
        "BCa: alphas %s -> %s (z0=%.4f, acc=%.4f)", alphas, new_alphas, z0, acc
    )
    return new_alphas


def segment_mean(cnarr: CopyNumArray, skip_low: bool = False) -> float64 | float:
    """Weighted average of bin log2 values."""
    if skip_low:
        cnarr = cnarr.drop_low_coverage()
    if len(cnarr) == 0:
        return np.nan
    if "weight" in cnarr and cnarr["weight"].any():
        return np.average(cnarr["log2"], weights=cnarr["weight"])  # type: ignore[no-any-return]
    return cnarr["log2"].mean()  # type: ignore[no-any-return]
