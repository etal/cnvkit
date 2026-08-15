"""Robust metrics to evaluate performance of copy number estimates."""

from __future__ import annotations

import hashlib
import logging
import warnings
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats

from . import descriptives

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

    from numpy import float64, ndarray
    from pandas.core.series import Series

    from cnvlib.cnary import CopyNumArray


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
            wt = weights[ser.index].to_numpy()
            vals = ser.to_numpy()
            # Drop bins with NaN weights to avoid poisoning np.average
            valid = ~np.isnan(wt)
            if valid.any():
                out_vals_lo[i], out_vals_hi[i] = func(vals[valid], wt[valid])
    return out_vals_lo, out_vals_hi


def confidence_interval_bootstrap(
    values: ndarray,
    weights: ndarray,
    alpha: float,
    bootstraps: int = 100,
    smoothed: bool | int = 10,
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

    # `bool` is a subclass of `int`, so the isinstance test cannot be folded into
    # the comparison: `smoothed=True` means "always smooth", not "threshold 1".
    use_smoothing = smoothed if isinstance(smoothed, bool) else k <= smoothed

    rng = _rng_from_values(values, weights)
    rand_indices = rng.integers(0, k, size=(bootstraps, k))
    samples: Iterable[tuple[ndarray, ndarray]] = (
        (np.take(values, idx), np.take(weights, idx)) for idx in rand_indices
    )
    if use_smoothing:
        samples = _smooth_samples_by_weight(k, samples, rng)
    # Recalculate segment means
    seg_means = (np.average(val, weights=wt) for val, wt in samples)
    bootstrap_dist = np.fromiter(seg_means, np.float64, bootstraps)
    alphas = np.array([alpha / 2, 1 - alpha / 2])
    if not use_smoothing:
        alphas = _bca_correct_alpha(values, weights, bootstrap_dist, alphas)
    ci = np.percentile(bootstrap_dist, list(100 * alphas))
    return ci  # type: ignore[no-any-return]


def _rng_from_values(values: ndarray, weights: ndarray) -> np.random.Generator:
    """Seed a generator from the segment's own values and weights.

    One fixed seed for every call gives every segment of the same bin count the
    same resample index matrix, so at the default 100 bootstraps the Monte Carlo
    error in ``ci_lo``/``ci_hi`` is a common shock across same-size segments
    rather than something that averages out along a profile: the bounds of two
    5-bin segments are correlated by construction. Mixing a digest of the
    segment's own numbers into the seed decorrelates them, and keeps both
    properties a fixed seed had:

    - repeated runs on the same input give identical bounds, since nothing
      outside the segment enters the seed;
    - an interval does not depend on the segment's position in the file or on
      what else the input contained, so ``segmetrics --ci`` on one chromosome
      reproduces that chromosome's rows of a whole-genome run. Keying on an
      iteration counter would have bought decorrelation at the cost of this.

    The digest covers the arrays as little-endian float64 bytes rather than as
    they arrive, pinning the byte order so the bounds do not depend on the
    machine. ``log2`` is float64 by construction (``CopyNumArray`` requires it),
    but ``weight`` keeps whatever dtype pandas inferred while reading the file,
    so an all-integral weight column arrives as int64 whose bytes are not those
    of the same numbers as float64 -- enough to break the subset invariance
    above, since which rows are integral depends on the slice. Both arrays have
    one entry per bin, so concatenating them cannot be read two ways. Python's
    ``hash`` is unusable here: it is salted by ``PYTHONHASHSEED``.
    """
    digest = hashlib.blake2b(digest_size=8)
    digest.update(np.asarray(values, dtype="<f8").tobytes())
    digest.update(np.asarray(weights, dtype="<f8").tobytes())
    seed = np.random.SeedSequence([0xA5EED, int.from_bytes(digest.digest(), "little")])
    return np.random.default_rng(seed)


def _smooth_samples_by_weight(
    k: int, samples: Iterable[tuple[ndarray, ndarray]], rng: np.random.Generator
) -> list[tuple[ndarray, ndarray]]:
    """Add Gaussian noise to each bootstrap replicate.

    The result is used to compute a "smoothed bootstrap," where the added noise
    ensures that for small samples (e.g. number of bins in the segment) the
    bootstrapped CI is close to the standard error of the mean, as it should be.
    Conceptually, sample from a KDE instead of the values themselves.

    This addresses the issue that small segments (#bins < #replicates) don't
    fully represent the underlying distribution, in particular the extreme
    values, so the CI is too narrow. Single-bin segments never arrive here:
    ``confidence_interval_bootstrap`` returns a zero-width interval for them
    before any resampling, so this widens the intervals of the multi-bin
    segments that ``smoothed`` selects.

    Standard deviation of the noise added to each bin comes from each bin's
    weight, which is an estimate of (1-variance).

    Parameters
    ----------
    k : int
        Number of bins in the segment.
    samples : iterable of (np.ndarray, np.ndarray)
        Bootstrap replicates as (value_sample, weight_sample).
    rng : np.random.Generator
        Seeded generator, shared with the caller's resampling draw so that the
        added noise is reproducible from run to run. Constructing a second
        generator here from the same seed would not make the two draws
        independent: identically seeded generators are the same stream.

    Returns
    -------
    `samples` with random N(0, pop_sd) added to each value, and
    weights unchanged.
    """
    # KDE bandwidth narrows for larger sample sizes, following Silverman's Rule
    # and Polansky 1995 but without their leading scale factor, which the
    # per-bin weight supplies instead.
    bw = k ** (-1 / 4)
    return [(v + (bw * np.sqrt(1 - w) * rng.standard_normal(k)), w) for v, w in samples]


def _jackknife_weighted_means(values: ndarray, weights: ndarray) -> ndarray:
    """Weighted mean of `values` with each element left out in turn.

    Taken from the totals rather than by re-averaging each leave-one-out subset:
    that is quadratic in the number of bins and dominated ``segmetrics --ci`` on
    WGS-scale profiles, where it accounted for 5.1 of the 5.3 seconds spent on a
    165k-bin fixture.

    Weights are non-negative, so the denominator vanishes only when every other
    bin has zero weight -- a case the subset average could not have handled
    either, since it raised on a zero weight sum. The resulting non-finite
    influence values propagate to a NaN alpha, which the caller's ``0 < a < 1``
    check rejects in favor of the uncorrected pair.
    """
    weighted = values * weights
    return (weighted.sum() - weighted) / (weights.sum() - weights)  # type: ignore[no-any-return]


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

    u = _jackknife_weighted_means(values, weights)
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
        # Drop bins with NaN weights to avoid poisoning np.average
        weights = cnarr["weight"].to_numpy()
        log2s = cnarr["log2"].to_numpy()
        valid = ~np.isnan(weights)
        if valid.any():
            return np.average(log2s[valid], weights=weights[valid])  # type: ignore[no-any-return]
    return cnarr["log2"].mean()  # type: ignore[no-any-return]
