"""Probe segmentation by convolving with the Haar wavelet.

The basic HaarSeg algorithm:

* Apply the undecimated discrete wavelet transform (UDWT) on the data, using the
  Haar wavelet.
* Select a set of detail subbands from the transform {LMIN, LMIN+1, ..., LMAX}.
* Find the local maxima of the selected detail subbands.
* Threshold the maxima of each subband separately, using an FDR thresholding
  procedure.
* Unify selected maxima from all the subbands to create a list of significant
  breakpoints in the data.
* Reconstruct the segmentation result from the list of significant breakpoints.

HaarSeg segmentation is based on detecting local maxima in the wavelet domain,
using Haar wavelet. The main algorithm parameter is breaksFdrQ, which controls
the sensitivity of the segmentation result. This function supports the optional
use of weights (also known as quality of measurments) and raw measurments. We
recommend using both extentions where possible, as it greatly improves the
segmentation result.

"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

if TYPE_CHECKING:
    from numpy import float64, ndarray

    from cnvlib.cnary import CopyNumArray
    from cnvlib.vary import VariantArray

# Lower bound for the wavelet-coefficient noise estimate, so a degenerate
# (all-flat / quantized) signal can't drive it to exactly 0 and break the FDR
# threshold's normal CDF. Far below any real log2/BAF step magnitude.
SIGMA_FLOOR = 1e-10

# Floor for per-bin BAF weights (read depths) so a window of all SNP-less bins
# can't divide by zero in the weighted Haar convolution. Negligible vs real
# read depths (~30-100x), so observed bins still dominate.
WEIGHT_FLOOR = 1e-6


def segment_haar(
    cnarr: CopyNumArray, fdr_q: float, variants: VariantArray | None = None
) -> CopyNumArray:
    """Do segmentation for CNVkit.

    Calculate copy number segmentation by HaarSeg
    (http://haarseg.r-forge.r-project.org/)

    Parameters
    ----------
    cnarr : CopyNumArray
        Binned, normalized copy ratios.
    fdr_q : float
        False discovery rate q-value.
    variants : VariantArray, optional
        Heterozygous SNP allele frequencies. When given, segment jointly on
        depth and B-allele frequency by taking the *union* of breakpoints from
        haar applied to log2 and to per-bin BAF -- a fast way to catch
        copy-neutral LOH on WGS (depth flat, BAF shifts). Without variants this
        is byte-identical to the depth-only segmentation.

    Returns
    -------
    CopyNumArray
        The CBS data table as a CNVkit object.
    """
    # Segment each chromosome individually
    # ENH - skip large gaps (segment chrom. arms separately)
    if variants is None:
        chrom_tables = [
            one_chrom(subprobes, fdr_q, chrom) for chrom, subprobes in cnarr.by_arm()
        ]
    else:
        chrom_tables = [
            one_chrom_baf(subprobes, fdr_q, chrom, variants)
            for chrom, subprobes in cnarr.by_arm()
        ]
    # ignore_index: by_arm splits a chromosome into arms, and each arm's segment
    # table starts its index at 0, so a plain concat would yield duplicate index
    # labels. A segmentation result must carry a unique index for downstream
    # label-based operations to behave (#1125).
    segarr = cnarr.as_dataframe(pd.concat(chrom_tables, ignore_index=True))
    segarr.sort_columns()
    return segarr


def _bin_weights(cnarr: CopyNumArray) -> ndarray | None:
    return cnarr["weight"].to_numpy() if "weight" in cnarr else None


def _haar_signal(cnarr: CopyNumArray) -> ndarray:
    """The raw per-bin log2 fed to haarSeg.

    HaarSeg is itself a denoising segmenter: its FDR threshold assumes the
    finest-scale (level-1) MAD reflects the true per-bin measurement noise.
    Pre-smoothing (e.g. Savitzky-Golay) deflates that noise floor while leaving
    coarse-scale autocorrelated structure intact, collapsing the FDR threshold
    so nearly every wiggle passes and ``fdr_q`` loses control -- ~8x
    over-segmentation into uniform ~7-bin pieces. Segmenting the
    raw log2 matches stock upstream HaarSeg and CNVkit's pre-2019 behavior; a
    prior scipy p-value bug (since fixed) was the real cause the smoothing was
    added to mask.
    """
    return cnarr["log2"].to_numpy()  # type: ignore[no-any-return]


def one_chrom(cnarr: CopyNumArray, fdr_q: float, chrom: str) -> pd.DataFrame:
    logging.debug("Segmenting %s", chrom)
    signal = _haar_signal(cnarr)
    weights = _bin_weights(cnarr)
    breakpoints = _haar_breakpoints(signal, fdr_q, weights)
    return _segments_from_breakpoints(cnarr, signal, weights, breakpoints, chrom)


def _haar_breakpoints(
    signal: ndarray,
    fdr_q: float,
    weights: ndarray | None = None,
    sigma_override: float | None = None,
) -> ndarray:
    """Bin-index breakpoints from haarSeg (i.e. ``results['start'][1:]``)."""
    results = haarSeg(signal, fdr_q, weights=weights, sigma_override=sigma_override)
    return results["start"][1:]


def _baf_signal_and_sigma(
    minor_counts: ndarray, depths: ndarray
) -> tuple[ndarray, ndarray, float]:
    """Build the BAF haar signal, per-bin weights, and noise floor.

    The signal is the pooled minor-allele fraction ``minor/depth`` per bin
    (already mirrored into [0, 0.5]); bins with no het SNP (depth 0) are
    neutral-filled at 0.5 and given weight 0. Weights are read depth, which is
    proportional to binomial precision (BAF variance ~= 0.25/depth).

    The noise floor is the MAD of level-1 Haar coefficients computed over
    *observed* bins only -- haarSeg's own ``peakSigmaEst`` is unweighted, so
    the neutral fills would otherwise deflate it and trip the SIGMA_FLOOR path.
    """
    observed = depths > 0
    baf_signal = np.where(observed, minor_counts / np.maximum(depths, 1.0), 0.5)
    baf_weights = depths.astype(np.float64)
    return baf_signal, baf_weights, _observed_baf_sigma(baf_signal, observed)


def _observed_baf_sigma(baf_signal: ndarray, observed: ndarray) -> float:
    """Noise floor from level-1 Haar coefficients of the *compacted* observed
    BAF values.

    Computing over the compacted observed values (not the full grid sampled at
    observed positions) means SNP-less neutral-fill bins neither deflate the
    estimate (flat 0.5 runs -> 0) nor inflate it (an observed LOH value vs its
    0.5-filled neighbors), which would otherwise raise the FDR threshold and
    suppress the LOH breakpoint we are trying to detect.
    """
    obs_vals = baf_signal[observed]
    if len(obs_vals) < 2:
        return SIGMA_FLOOR
    level1 = HaarConv(obs_vals, None, 1)
    return max(float(pd.Series(level1).abs().median() * 1.4826), SIGMA_FLOOR)


def _baf_signal_from_frequencies(
    baf_series: pd.Series,
) -> tuple[ndarray, ndarray, float]:
    """Fallback BAF signal/weights/noise when the VCF lacks allele counts.

    Uses the per-bin mirrored BAF (NaN where no SNP) with binary observed/not
    weights, and the same observed-only noise estimate as the count-based path.
    """
    vals = baf_series.to_numpy(dtype=np.float64)
    observed = ~np.isnan(vals)
    baf_signal = np.where(observed, vals, 0.5)
    baf_weights = observed.astype(np.float64)
    return baf_signal, baf_weights, _observed_baf_sigma(baf_signal, observed)


def one_chrom_baf(
    cnarr: CopyNumArray, fdr_q: float, chrom: str, variants: VariantArray
) -> pd.DataFrame:
    """Joint depth+BAF haar segmentation of one chromosome arm.

    Union the breakpoints from haar on log2 (depth) and haar on per-bin BAF,
    then rebuild segments. Detects copy-neutral LOH (BAF shift, flat depth).
    """
    logging.debug("Segmenting (depth+BAF) %s", chrom)
    signal = _haar_signal(cnarr)
    weights = _bin_weights(cnarr)
    log2_bps = _haar_breakpoints(signal, fdr_q, weights)

    # Build the BAF signal: prefer count-aware (depth-weighted) over frequencies,
    # and skip BAF entirely if the VCF carries neither counts nor frequencies.
    counts = variants.baf_counts_by_ranges(cnarr)
    if counts is not None:
        minor_counts, depths = (c.to_numpy(dtype=np.float64) for c in counts)
        baf_signal, baf_weights, baf_sigma = _baf_signal_and_sigma(minor_counts, depths)
    elif "alt_freq" in variants:
        baf_signal, baf_weights, baf_sigma = _baf_signal_from_frequencies(
            variants.baf_by_ranges(cnarr)
        )
    else:
        baf_weights = None  # no usable allele data -> depth-only below

    if baf_weights is not None and np.count_nonzero(baf_weights) >= 2:
        baf_bps = _haar_breakpoints(
            baf_signal,
            fdr_q,
            weights=np.maximum(baf_weights, WEIGHT_FLOOR),
            sigma_override=baf_sigma,
        )
        breakpoints = UnifyLevels(log2_bps, baf_bps, 1)
    else:
        # Too few informative bins on this arm -> depth-only (no BAF breaks)
        breakpoints = log2_bps
    return _segments_from_breakpoints(cnarr, signal, weights, breakpoints, chrom)


def _segments_from_breakpoints(
    cnarr: CopyNumArray,
    signal: ndarray,
    weights: ndarray | None,
    breakpoints: ndarray,
    chrom: str,
) -> pd.DataFrame:
    """Build a segment table from breakpoint bin indices (shared by haar paths).

    `signal` is the raw per-bin log2 (from `_haar_signal`) used to find the
    breakpoints; reuse the same array here so the segment means are consistent
    with the breakpoints.
    """
    n = len(cnarr)
    seg_starts = np.insert(breakpoints, 0, 0)
    seg_ends = np.append(breakpoints, n)
    log2_means = SegmentByPeaks(signal, breakpoints, weights)
    return pd.DataFrame(
        {
            "chromosome": chrom,
            "start": cnarr["start"].to_numpy().take(seg_starts),
            "end": cnarr["end"].to_numpy().take(seg_ends - 1),
            "log2": log2_means[seg_starts],
            "gene": "-",
            "probes": seg_ends - seg_starts,
        }
    )


def haarSeg(
    signal: ndarray,
    breaksFdrQ: float,
    weights: ndarray | None = None,
    raw_signal: ndarray | None = None,
    haarStartLevel: int = 1,
    haarEndLevel: int = 5,
    sigma_override: float | None = None,
) -> dict[str, ndarray]:
    r"""Perform segmentation according to the HaarSeg algorithm.

    Parameters
    ----------
    signal : array
        A 1D array of log-ratio values, sorted according to their genomic
        location.
    weights : array
        Weight matrix, corresponding to quality of measurement, with values
        :math:`1/(\sigma^2)`. Must have the same size as signal.
    raw_signal : array
        The minimum between the raw test-sample and control-sample coverages
        (before applying log ratio, but after any background reduction and/or
        normalization). These raw red / green measurments are used to detect
        low-value probes, which are more sensitive to noise.
        Used for the non-stationary variance compensation.
        Must have the same size as signal.
    breaksFdrQ : float
        The FDR q parameter. This value should lie between 0 and 0.5. The
        smaller this value is, the less sensitive the segmentation result will
        be.
        For example, we will detect fewer segmentation breaks when using Q =
        1e-4, compared to when using Q = 1e-3.
        Common used values are 1e-2, 1e-3, 1e-4.
    haarStartLevel : int
        The detail subband from which we start to detect peaks. The higher this
        value is, the less sensitive we are to short segments. The default is
        value is 1, corresponding to segments of 2 probes.
    haarEndLevel : int
        The detail subband until which we use to detect peaks. The higher this
        value is, the more sensitive we are to large trends in the data. This
        value DOES NOT indicate the largest possible segment that can be
        detected.  The default is value is 5, corresponding to step of 32 probes
        in each direction.

    Returns
    -------
    dict

    Source: haarSeg.R
    """

    def med_abs_diff(diff_vals):
        """Median absolute deviation, with deviations given."""
        if len(diff_vals) == 0:
            return 0.0
        return diff_vals.abs().median() * 1.4826

    diff_signal = pd.Series(HaarConv(signal, None, 1))
    if raw_signal:  # type: ignore[unreachable]
        # Non-stationary variance empirical threshold set to 50
        NSV_TH = 50
        varMask = raw_signal < NSV_TH
        pulseSize = 2
        diffMask = PulseConv(varMask, pulseSize) >= 0.5
        peakSigmaEst = med_abs_diff(diff_signal[~diffMask])
        noisySigmaEst = med_abs_diff(diff_signal[diffMask])
    else:
        peakSigmaEst = med_abs_diff(diff_signal)

    # Floor the noise estimate: for mostly-flat / quantized input (>=50% of
    # adjacent bins bit-identical) the median absolute coefficient is exactly 0,
    # which would make FDRThres' norm.cdf(scale=0) return NaN p-values and
    # silently drop real breakpoints. A tiny positive floor keeps genuine steps
    # (magnitude >> SIGMA_FLOOR) significant without admitting numerical noise.
    # Real continuous log2/BAF data has noise, so peakSigmaEst >> the floor and
    # this is a no-op there.
    peakSigmaEst = max(peakSigmaEst, SIGMA_FLOOR)
    if sigma_override is not None:
        # Caller supplies a noise estimate computed over informative bins only
        # (the BAF pass: peakSigmaEst above is contaminated by neutral-filled,
        # SNP-less bins because it is computed unweighted). Default None keeps
        # the depth pass byte-identical.
        peakSigmaEst = max(sigma_override, SIGMA_FLOOR)

    breakpoints = np.array([], dtype=np.int_)
    for level in range(haarStartLevel, haarEndLevel + 1):
        stepHalfSize = 2**level
        convRes = HaarConv(signal, weights, stepHalfSize)
        peakLoc = FindLocalPeaks(convRes)
        logging.debug("Found %d peaks at level %d", len(peakLoc), level)

        if raw_signal:  # type: ignore[unreachable]
            pulseSize = 2 * stepHalfSize
            convMask = PulseConv(varMask, pulseSize) >= 0.5
            sigmaEst = (1 - convMask) * peakSigmaEst + convMask * noisySigmaEst
            convRes /= sigmaEst
            peakSigmaEst = 1.0

        T = FDRThres(convRes[peakLoc], breaksFdrQ, peakSigmaEst)
        # Keep only the peak values where the signal amplitude is large enough.
        addonPeaks = np.extract(np.abs(convRes.take(peakLoc)) >= T, peakLoc)
        breakpoints = UnifyLevels(breakpoints, addonPeaks, 2 ** (level - 1))

    logging.debug("Found %d breakpoints: %s", len(breakpoints), breakpoints)

    # Translate breakpoints to segments
    segs = SegmentByPeaks(signal, breakpoints, weights)
    segSt = np.insert(breakpoints, 0, 0)
    segEd = np.append(breakpoints, len(signal))
    return {
        "start": segSt,
        "end": segEd - 1,
        "size": segEd - segSt,
        "mean": segs[segSt],
    }


def FDRThres(x: ndarray, q: float, stdev: float64) -> int | float64:
    """False discovery rate (FDR) threshold."""
    M = len(x)
    if M < 2:
        return 0
    if stdev <= 0:
        # Zero/negative noise: norm.cdf(scale<=0) is NaN. Against a zero noise
        # floor every nonzero peak is infinitely significant, so admit them all
        # (threshold 0); the caller floors stdev (SIGMA_FLOOR) so this is a guard.
        return np.float64(0.0)

    m = np.arange(1, M + 1) / M
    x_sorted = np.sort(np.abs(x))[::-1]
    # Two-tailed p-value: P(|X| > x) where X ~ N(0, stdev)
    p = 2 * (1 - stats.norm.cdf(x_sorted, loc=0, scale=stdev))
    # Get the largest index for which p <= m*q
    indices = np.nonzero(p <= m * q)[0]
    if len(indices):
        T = x_sorted[indices[-1]]
    else:
        logging.debug(
            "No passing p-values: min p=%.4g, min m=%.4g, q=%s", p[0], m[0], q
        )
        # No peak is significant -> set the threshold just above the largest so
        # the '>=' test admits nothing. Use nextafter, not '+1e-16': near
        # magnitude 1 that addend is below the float64 ULP (~2.22e-16) and is a
        # no-op, which let the largest peak slip through.
        T = np.nextafter(x_sorted[0], np.inf)
    return T  # type: ignore[no-any-return]


def SegmentByPeaks(
    data: ndarray, peaks: ndarray, weights: ndarray | None = None
) -> ndarray:
    """Average the values of the probes within each segment.

    Parameters
    ----------
    data : array
        the probe array values
    peaks : array
        Positions of copy number breakpoints in the original array

    Source: SegmentByPeaks.R
    """
    if len(peaks) == 0:
        # Single segment spanning all data
        if weights is not None and np.nansum(weights) > 0:
            valid = ~np.isnan(weights)
            mean_val = np.average(data[valid], weights=weights[valid])
        else:
            mean_val = np.mean(data)
        return np.full_like(data, mean_val)

    # Segment boundaries: [0, peaks..., len(data)]
    bounds = np.concatenate([[0], peaks, [len(data)]])
    seg_lengths = np.diff(bounds)
    n_segs = len(seg_lengths)

    if weights is None:
        # Vectorized unweighted means using reduceat
        seg_sums = np.add.reduceat(data, bounds[:-1])
        seg_means = seg_sums / seg_lengths
    else:
        # Weighted means require per-segment computation
        seg_means = np.empty(n_segs, dtype=np.float64)
        for i in range(n_segs):
            start, end = bounds[i], bounds[i + 1]
            seg_data = data[start:end]
            w = weights[start:end]
            if np.nansum(w) > 0:
                valid = ~np.isnan(w)
                seg_means[i] = np.average(seg_data[valid], weights=w[valid])
            else:
                seg_means[i] = np.mean(seg_data)

    # Expand segment means to full array using repeat
    return np.repeat(seg_means, seg_lengths)


# ---- from HaarSeg C code -- the core ----


# --- HaarSeg.h


def HaarConv(
    signal: ndarray,  # const double * signal,
    weight: ndarray | None,  # const double * weight,
    stepHalfSize: int,  # int stepHalfSize,
) -> ndarray:
    """Convolve haar wavelet function with a signal, applying circular padding.

    Parameters
    ----------
    signal : const array of floats
    weight : const array of floats (optional)
    stepHalfSize : int

    Returns
    -------
    array
        Of floats, representing the convolved signal.

    Source: HaarSeg.c
    """
    signalSize = len(signal)
    if stepHalfSize > signalSize:
        # Signal too short for this wavelet scale; return zeros to skip this level
        logging.warning(
            "Wavelet step size (%d) exceeds signal length (%d); "
            "skipping this decomposition level",
            stepHalfSize,
            signalSize,
        )
        return np.zeros(signalSize, dtype=np.float64)

    if weight is None:
        # Vectorized unweighted case
        return _haar_conv_unweighted(signal, stepHalfSize)
    else:
        # Weighted case (not vectorized due to running sums with dependencies)
        return _haar_conv_weighted(signal, weight, stepHalfSize)


def _haar_conv_unweighted(signal: ndarray, stepHalfSize: int) -> ndarray:
    """Vectorized Haar convolution for unweighted signals."""
    signalSize = len(signal)
    k = np.arange(1, signalSize)

    # Compute highEnd indices with circular padding
    highEnd = k + stepHalfSize - 1
    over_mask = highEnd >= signalSize
    highEnd[over_mask] = signalSize - 1 - (highEnd[over_mask] - signalSize)

    # Compute lowEnd indices with circular padding
    lowEnd = k - stepHalfSize - 1
    under_mask = lowEnd < 0
    lowEnd[under_mask] = -lowEnd[under_mask] - 1

    # Compute increments and cumulative sum
    increments = signal[highEnd] + signal[lowEnd] - 2 * signal[k - 1]
    result = np.zeros(signalSize, dtype=np.float64)
    result[1:] = np.cumsum(increments)

    # Normalize
    stepNorm = math.sqrt(2.0 * stepHalfSize)
    result[1:] /= stepNorm

    return result


def _haar_conv_weighted(signal: ndarray, weight: ndarray, stepHalfSize: int) -> ndarray:
    """Haar convolution for weighted signals (sequential due to running sums)."""
    signalSize = len(signal)
    result = np.zeros(signalSize, dtype=np.float64)

    # Init weight sums
    highWeightSum = weight[:stepHalfSize].sum()
    highNonNormed = (weight[:stepHalfSize] * signal[:stepHalfSize]).sum()
    # Circular padding
    lowWeightSum = highWeightSum
    lowNonNormed = -highNonNormed

    norm_factor = math.sqrt(stepHalfSize / 2)

    for k in range(1, signalSize):
        highEnd = k + stepHalfSize - 1
        if highEnd >= signalSize:
            highEnd = signalSize - 1 - (highEnd - signalSize)
        lowEnd = k - stepHalfSize - 1
        if lowEnd < 0:
            lowEnd = -lowEnd - 1

        lowNonNormed += signal[lowEnd] * weight[lowEnd] - signal[k - 1] * weight[k - 1]
        highNonNormed += (
            signal[highEnd] * weight[highEnd] - signal[k - 1] * weight[k - 1]
        )
        lowWeightSum += weight[k - 1] - weight[lowEnd]
        highWeightSum += weight[highEnd] - weight[k - 1]
        result[k] = norm_factor * (
            lowNonNormed / lowWeightSum + highNonNormed / highWeightSum
        )

    return result


def FindLocalPeaks(
    signal: ndarray,  # const double * signal,
    # peakLoc, #int * peakLoc
) -> ndarray:
    """Find local maxima on positive values, local minima on negative values.

    First and last index are never considered extramum.

    Parameters
    ----------
    signal : const array of floats

    Returns
    -------
    peakLoc : array of ints
        Locations of extrema in `signal`

    Source: HaarSeg.c
    """
    # use numpy.diff to simplify? argmax, argmin?
    maxSuspect = minSuspect = None
    peakLoc = []
    for k in range(1, len(signal) - 1):
        sig_prev, sig_curr, sig_next = signal[k - 1 : k + 2]
        if sig_curr > 0:
            # Look for local maxima
            if (sig_curr > sig_prev) and (sig_curr > sig_next):
                peakLoc.append(k)
            elif (sig_curr > sig_prev) and (sig_curr == sig_next):
                maxSuspect = k
            elif (sig_curr == sig_prev) and (sig_curr > sig_next):
                # Take the first in a series of equal values
                if maxSuspect is not None:
                    peakLoc.append(maxSuspect)
                    maxSuspect = None
            elif (sig_curr == sig_prev) and (sig_curr < sig_next):
                maxSuspect = None

        elif sig_curr < 0:
            # Look for local maxima
            if (sig_curr < sig_prev) and (sig_curr < sig_next):
                peakLoc.append(k)
            elif (sig_curr < sig_prev) and (sig_curr == sig_next):
                minSuspect = k
            elif (sig_curr == sig_prev) and (sig_curr < sig_next):
                if minSuspect is not None:
                    peakLoc.append(minSuspect)
                    minSuspect = None
            elif (sig_curr == sig_prev) and (sig_curr > sig_next):
                minSuspect = None

    return np.array(peakLoc, dtype=np.int_)


def UnifyLevels(
    baseLevel: ndarray,  # const int * baseLevel,
    addonLevel: ndarray,  # const int * addonLevel,
    windowSize: int,  # int windowSize,
) -> ndarray:
    """Unify several decomposition levels.

    Merge the two lists of breakpoints, but drop addonLevel values that are too
    close to baseLevel values.

    Parameters
    ----------
    baseLevel : const array of ints
    addonLevel : const array of ints
    windowSize : int

    Returns
    -------
    joinedLevel : array of ints

    Source: HaarSeg.c
    """
    if not len(addonLevel):
        return baseLevel
    if not len(baseLevel):
        return addonLevel.copy()  # type: ignore[no-any-return]

    # Use searchsorted to find nearest base elements for each addon
    # Since baseLevel is sorted, the closest base is either at insert_pos or insert_pos-1
    insert_pos = np.searchsorted(baseLevel, addonLevel)
    n_base = len(baseLevel)

    # Check distance to left neighbor (where it exists)
    has_left = insert_pos > 0
    left_idx = np.clip(insert_pos - 1, 0, n_base - 1)
    left_too_close = has_left & ((addonLevel - baseLevel[left_idx]) <= windowSize)

    # Check distance to right neighbor (where it exists)
    has_right = insert_pos < n_base
    right_idx = np.clip(insert_pos, 0, n_base - 1)
    right_too_close = has_right & ((baseLevel[right_idx] - addonLevel) <= windowSize)

    # Keep addon elements that are not too close to any base element
    keep_mask = ~(left_too_close | right_too_close)

    # Merge kept addons with base and sort
    joined = np.concatenate([baseLevel, addonLevel[keep_mask]])
    joined.sort()
    return joined  # type: ignore[no-any-return]


def PulseConv(
    signal: ndarray,
    pulseSize: int,
) -> ndarray:
    """Convolve a pulse function with a signal, applying circular padding.

    Used for non-stationary variance compensation.

    Parameters
    ----------
    signal : array of floats
    pulseSize : int

    Returns
    -------
    array of floats

    Source: HaarSeg.c
    """
    signalSize = len(signal)
    if pulseSize > signalSize:
        # Signal too short for this pulse size; return zeros to skip
        logging.warning(
            "Pulse size (%d) exceeds signal length (%d); skipping convolution",
            pulseSize,
            signalSize,
        )
        return np.zeros(signalSize, dtype=np.float64)

    pulseHeight = 1.0 / pulseSize

    # Compute initial value with circular padding
    result_0 = signal[: (pulseSize + 1) // 2].sum() + signal[: pulseSize // 2].sum()
    result_0 *= pulseHeight

    # Vectorized main computation
    k = np.arange(pulseSize // 2, signalSize + (pulseSize // 2) - 1)

    # Compute tail indices with circular padding
    tail = k - pulseSize
    tail_neg = tail < 0
    tail[tail_neg] = -tail[tail_neg] - 1

    # Compute head indices with circular padding
    head = k.copy()
    head_over = head >= signalSize
    head[head_over] = signalSize - 1 - (head[head_over] - signalSize)

    # Compute increments and cumulative sum
    increments = (signal[head] - signal[tail]) * pulseHeight
    result = np.empty(signalSize, dtype=np.float64)
    result[0] = result_0
    result[1:] = result_0 + np.cumsum(increments)

    return result


# XXX Apply afterward to the segmentation result? (not currently used)


def AdjustBreaks(
    signal,  # const double * signal,
    peakLoc,  # const int * peakLoc,
):
    """Improve localization of breaks. Suboptimal, but linear-complexity.

    We try to move each break 1 sample left/right, choosing the offset which
    leads to minimum data error.

    Parameters
    ----------
    signal: const array of floats
    peakLoc: const array of ints

    Source: HaarSeg.c
    """
    newPeakLoc = peakLoc.copy()
    for k, npl_k in enumerate(newPeakLoc):
        # Calculating width of segments around the breakpoint
        n1 = npl_k if k == 0 else npl_k - newPeakLoc[k - 1]
        n2 = (len(signal) if k + 1 == len(newPeakLoc) else newPeakLoc[k + 1]) - npl_k

        # Find the best offset for current breakpoint, trying only 1 sample
        # offset
        bestScore = float("Inf")  # Smaller is better
        bestOffset = 0
        for p in (-1, 0, 1):
            # Pointless to try to remove single-sample segments
            if (n1 == 1 and p == -1) or (n2 == 1 and p == 1):
                continue

            signal_n1_to_p = signal[npl_k - n1 : npl_k + p]
            s1 = signal_n1_to_p.sum() / (n1 + p)
            ss1 = ((signal_n1_to_p - s1) ** 2).sum()

            signal_p_to_n2 = signal[npl_k + p : npl_k + n2]
            s2 = signal_p_to_n2.sum() / (n2 - p)
            ss2 = ((signal_p_to_n2 - s2) ** 2).sum()

            score = ss1 + ss2
            if score < bestScore:
                bestScore = score
                bestOffset = p

        if bestOffset != 0:
            newPeakLoc[k] += bestOffset

    return newPeakLoc


# Testing
