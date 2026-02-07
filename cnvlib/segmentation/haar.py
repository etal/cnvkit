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
    from cnvlib.cnary import CopyNumArray
    from numpy import float64, ndarray


def segment_haar(cnarr: CopyNumArray, fdr_q: float) -> CopyNumArray:
    """Do segmentation for CNVkit.

    Calculate copy number segmentation by HaarSeg
    (http://haarseg.r-forge.r-project.org/)

    Parameters
    ----------
    cnarr : CopyNumArray
        Binned, normalized copy ratios.
    fdr_q : float
        False discovery rate q-value.

    Returns
    -------
    CopyNumArray
        The CBS data table as a CNVkit object.
    """
    # Segment each chromosome individually
    # ENH - skip large gaps (segment chrom. arms separately)
    chrom_tables = [
        one_chrom(subprobes, fdr_q, chrom) for chrom, subprobes in cnarr.by_arm()
    ]
    segarr = cnarr.as_dataframe(pd.concat(chrom_tables))
    segarr.sort_columns()
    return segarr


def one_chrom(cnarr: CopyNumArray, fdr_q: float, chrom: str) -> pd.DataFrame:
    logging.debug("Segmenting %s", chrom)
    results = haarSeg(
        cnarr.smooth_log2(),
        fdr_q,
        weights=(cnarr["weight"].to_numpy() if "weight" in cnarr else None),
    )
    table = pd.DataFrame(
        {
            "chromosome": chrom,
            "start": cnarr["start"].to_numpy().take(results["start"]),
            "end": cnarr["end"].to_numpy().take(results["end"]),
            "log2": results["mean"],
            "gene": "-",
            "probes": results["size"],
        }
    )
    return table


def variants_in_segment(varr, segment, fdr_q):
    """Segment a single genomic interval based on B-allele frequencies.

    Applies HaarSeg segmentation to variant allele frequencies within a segment
    to detect sub-clonal changes or allelic imbalances. This is used for
    allele-specific copy number analysis.

    Parameters
    ----------
    varr : VariantArray
        Variant allele frequency data (from VCF) within the segment region.
        Contains SNV positions and their B-allele frequencies.
    segment : Row or CopyNumArray
        A single segment from copy number segmentation results.
        Must have 'chromosome', 'start', 'end', 'gene', 'log2', 'probes' fields.
    fdr_q : float
        False discovery rate q-value for HaarSeg breakpoint detection.
        Typical values: 0.01, 0.001, 0.0001. Lower = less sensitive.

    Returns
    -------
    pandas.DataFrame
        Segmentation results as a table with columns:
        - chromosome, start, end: Genomic coordinates
        - gene: Gene name from parent segment
        - log2: Copy ratio from parent segment (not re-calculated)
        - probes: Number of variants in each sub-segment

        If no sub-segmentation is detected (≤1 breakpoint), returns the
        original segment as a single-row DataFrame.

    Notes
    -----
    The function:

    1. Transforms B-allele frequencies to mirrored values (0.5-1.0 range)
       with tumor purity boosting
    2. Applies HaarSeg to detect breakpoints in allele frequencies
    3. If multiple segments found, places breakpoints midway between SNVs
    4. If no segmentation needed, returns original segment unchanged

    The log2 copy ratio is inherited from the parent segment and not
    recalculated, as this function focuses on allelic changes.

    This is primarily used internally by the `call` command when VCF data
    is provided to refine segment boundaries based on heterozygous SNPs.

    See Also
    --------
    haarSeg : The underlying HaarSeg segmentation algorithm
    VariantArray.mirrored_baf : Transforms BAF values for segmentation
    """
    if len(varr):
        values = varr.mirrored_baf(above_half=True, tumor_boost=True)
        results = haarSeg(values, fdr_q, weights=None)  # ENH weight by sqrt(DP)
    else:
        values = pd.Series()
        results = None
    if results is not None and len(results["start"]) > 1:
        logging.info(
            "Segmented on allele freqs in %s:%d-%d",
            segment.chromosome,
            segment.start,
            segment.end,
        )
        # Ensure breakpoint locations make sense
        # - Keep original segment start, end positions
        # - Place breakpoints midway between SNVs, I guess?
        # NB: 'results' are indices, i.e. enumerated bins
        gap_rights = varr["start"].to_numpy().take(results["start"][1:])
        gap_lefts = varr["end"].to_numpy().take(results["end"][:-1])
        mid_breakpoints = (gap_lefts + gap_rights) // 2
        starts = np.concatenate([[segment.start], mid_breakpoints])
        ends = np.concatenate([mid_breakpoints, [segment.end]])
        table = pd.DataFrame(
            {
                "chromosome": segment.chromosome,
                "start": starts,
                "end": ends,
                # 'baf': results['mean'],
                "gene": segment.gene,  # '-'
                "log2": segment.log2,
                "probes": results["size"],
                # 'weight': (segment.weight * results['size']
                #            / (segment.end - segment.start)),
            }
        )
    else:
        table = pd.DataFrame(
            {
                "chromosome": segment.chromosome,
                "start": segment.start,
                "end": segment.end,
                # 'baf': values.median(),
                "gene": segment.gene,  # '-',
                "log2": segment.log2,
                "probes": segment.probes,
                # 'weight': segment.weight,
            },
            index=[0],
        )

    return table


# ---- from HaarSeg R code -- the API ----


def haarSeg(
    signal: ndarray,
    breaksFdrQ: float,
    weights: ndarray | None = None,
    raw_signal: ndarray | None = None,
    haarStartLevel: int = 1,
    haarEndLevel: int = 5,
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
        T = x_sorted[0] + 1e-16  # ~= 2^-52, like MATLAB "eps"
    return T


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
        if weights is not None and weights.sum() > 0:
            mean_val = np.average(data, weights=weights)
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
            w = weights[start:end]
            if w.sum() > 0:
                seg_means[i] = np.average(data[start:end], weights=w)
            else:
                seg_means[i] = np.mean(data[start:end])

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
        return addonLevel.copy()

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
    return joined


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


def table2coords(seg_table):
    """Return x, y arrays for plotting."""
    x = []
    y = []
    for start, size, val in seg_table:
        x.append(start)
        x.append(start + size)
        y.append(val)
        y.append(val)
    return x, y


if __name__ == "__main__":
    real_data = np.concatenate(
        (np.zeros(800), np.ones(200), np.zeros(800), 0.8 * np.ones(200), np.zeros(800))
    )
    # rng = np.random.default_rng(0x5EED)
    rng = np.random.default_rng()
    noisy_data = real_data + rng.standard_normal(len(real_data)) * 0.2

    # # Run using default parameters
    seg_table = haarSeg(noisy_data, 0.005)

    logging.info("%s", seg_table)

    from matplotlib import pyplot

    indices = np.arange(len(noisy_data))
    pyplot.scatter(indices, noisy_data, alpha=0.2, color="gray")
    x, y = table2coords(seg_table)
    pyplot.plot(x, y, color="r", marker="x", lw=2, snap=False)
    pyplot.show()

    # # The complete segmented signal
    # lines(seg.data$Segmented, col="red", lwd=3)
