#!/usr/bin/env python

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
from __future__ import absolute_import, division, print_function
from builtins import range, zip

import logging
import math

import numpy as np
import pandas as pd
from scipy import stats


def segment_haar(cnarr, fdr_q):
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
    chrom_tables = [one_chrom(subprobes, fdr_q, chrom)
                    for chrom, subprobes in cnarr.by_chromosome()]
    segarr = cnarr.as_dataframe(pd.concat(chrom_tables))
    segarr.sort_columns()
    return segarr


def one_chrom(cnarr, fdr_q, chrom):
    logging.debug("Segmenting %s", chrom)
    segtable = haarSeg(np.asarray(cnarr['log2']),
                       fdr_q,
                       W=(np.asarray(cnarr['weight'])
                          if 'weight' in cnarr
                          else None))
    table = pd.DataFrame({
        'chromosome': chrom,
        'start': np.asarray(cnarr['start']).take(segtable['start']),
        'end': np.asarray(cnarr['end']).take(segtable['end']),
        'log2': segtable['mean'],
        'gene': '-',
        'probes': segtable['size'],
    })
    return table


def variants_in_segment(varr, segment, fdr_q):
    if len(varr):
        values = varr.mirrored_baf(above_half=True, tumor_boost=True)
        segtable = haarSeg(values,
                           fdr_q,
                           W=None)  # weight by sqrt(DP)?
    else:
        values = pd.Series()
        segtable = None
    if segtable is not None and len(segtable['start']) > 1:
        logging.info("Segmented on allele freqs in %s:%d-%d",
                     segment.chromosome, segment.start, segment.end)
        # Ensure breakpoint locations make sense
        # - Keep original segment start, end positions
        # - Place breakpoints midway between SNVs, I guess?
        gap_rights = np.asarray(varr['start']).take(segtable['start'][1:])
        gap_lefts = np.asarray(varr['end']).take(segtable['end'][:-1])
        mid_breakpoints = [(left + right) // 2
                           for left, right in zip(gap_lefts, gap_rights)]
        starts = np.concatenate([[segment.start], mid_breakpoints])
        ends = np.concatenate([mid_breakpoints, [segment.end]])
        table = pd.DataFrame({
            'chromosome': segment.chromosome,
            'start': starts,
            'end': ends,
            # 'baf': segtable['mean'],
            'gene': segment.gene, # '-'
            'log2': segment.log2,
            'probes': segtable['size'],
            # 'weight': (segment.weight * segtable['size']
            #            / (segment.end - segment.start)),
        })
    else:
        table = pd.DataFrame({
            'chromosome': segment.chromosome,
            'start': segment.start,
            'end': segment.end,
            # 'baf': values.median(),
            'gene': segment.gene, #'-',
            'log2': segment.log2,
            'probes': segment.probes,
            # 'weight': segment.weight,
        }, index=[0])

    return table



# ---- from HaarSeg R code -- the API ----

def haarSeg(I, breaksFdrQ,
            W=None,
            rawI=None,
            haarStartLevel=1,
            haarEndLevel=5):
    r"""Perform segmentation according to the HaarSeg algorithm.

    Parameters
    ----------
    I : array
        A 1D array of log-ratio values, sorted according to their genomic
        location.
    W : array
        Weight matrix, corresponding to quality of measurement, with values
        :math:`1/(\sigma^2)`. Must have the same size as I.
    rawI : array
        The minimum between the raw test-sample and control-sample coverages
        (before applying log ratio, but after any background reduction and/or
        normalization). These raw red / green measurments are used to detect
        low-value probes, which are more sensitive to noise.
        Used for the non-stationary variance compensation.
        Must have the same size as I.
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
    tuple
        Two elements:
        1. Segments result table, a list of tuples:
            (segment start index, segment size, segment value)
        2. Segmented signal array (same size as I)

    Source: haarSeg.R
    """
    def med_abs_diff(diff_vals):
        """Median absolute deviation, with deviations given."""
        if len(diff_vals) == 0:
            return 0.
        return diff_vals.abs().median() * 1.4826

    diffI = pd.Series(HaarConv(I, None, 1))
    if rawI:
        # Non-stationary variance empirical threshold set to 50
        NSV_TH = 50
        varMask = (rawI < NSV_TH)
        pulseSize = 2
        diffMask = (PulseConv(varMask, pulseSize) >= .5)
        peakSigmaEst = med_abs_diff(diffI[~diffMask])
        noisySigmaEst = med_abs_diff(diffI[diffMask])
    else:
        peakSigmaEst = med_abs_diff(diffI)

    breakpoints = np.array([], dtype=np.int_)
    for level in range(haarStartLevel, haarEndLevel+1):
        stepHalfSize = 2 ** level
        convRes = HaarConv(I, W, stepHalfSize)
        peakLoc = FindLocalPeaks(convRes)
        logging.debug("Found %d peaks at level %d", len(peakLoc), level)

        if rawI:
            pulseSize = 2 * stepHalfSize
            convMask = (PulseConv(varMask, pulseSize) >= .5)
            sigmaEst = (1 - convMask) * peakSigmaEst + convMask * noisySigmaEst
            convRes /= sigmaEst
            peakSigmaEst = 1.

        T = FDRThres(convRes[peakLoc], breaksFdrQ, peakSigmaEst)
        # Keep only the peak values where the signal amplitude is large enough.
        addonPeaks = np.extract(np.abs(convRes.take(peakLoc)) >= T, peakLoc)
        breakpoints = UnifyLevels(breakpoints, addonPeaks, 2 ** (level - 1))

    logging.debug("Found %d breakpoints: %s", len(breakpoints), breakpoints)

    # Translate breakpoints to segments
    segs = SegmentByPeaks(I, breakpoints, W)
    segSt = np.insert(breakpoints, 0, 0)
    segEd = np.append(breakpoints, len(I))
    return {'start': segSt,
            'end': segEd - 1,
            'size': segEd - segSt,
            'mean': segs[segSt]}


def FDRThres(x, q, stdev):
    """False discovery rate (FDR) threshold."""
    M = len(x)
    if M < 2:
        return 0

    m = np.arange(1, M+1) / M
    x_sorted = np.sort(np.abs(x))[::-1]
    p = 2 * (1 - stats.norm.cdf(x_sorted, stdev))  # like R "pnorm"
    # Get the largest index for which p <= m*q
    indices = np.nonzero(p <= m * q)[0]
    if len(indices):
        T = x_sorted[indices[-1]]
    else:
        logging.debug("No passing p-values: min p=%.4g, min m=%.4g, q=%s",
                      p[0], m[0], q)
        T = x_sorted[0] + 1e-16  # ~= 2^-52, like MATLAB "eps"
    return T


def SegmentByPeaks(data, peaks, weights=None):
    """Average the values of the probes within each segment.

    Parameters
    ----------
    data : array
        the probe array values
    peaks : array
        Positions of copy number breakpoints in the original array

    Source: SegmentByPeaks.R
    """
    segs = np.zeros_like(data)
    for seg_start, seg_end in zip(np.insert(peaks, 0, 0),
                                  np.append(peaks, len(data))):
        if weights is not None and weights[seg_start:seg_end].sum() > 0:
            # Weighted mean of individual probe values
            val = np.average(data[seg_start:seg_end],
                             weights=weights[seg_start:seg_end])
        else:
            # Unweighted mean of individual probe values
            val = np.mean(data[seg_start:seg_end])
        segs[seg_start:seg_end] = val
    return segs



# ---- from HaarSeg C code -- the core ----

# --- HaarSeg.h
def HaarConv(signal, #const double * signal,
             weight, #const double * weight,
             stepHalfSize, #int stepHalfSize,
            ):
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
        # XXX TODO handle this endcase
        # raise ValueError("stepHalfSize (%s) > signalSize (%s)"
        #                  % (stepHalfSize, signalSize))
        logging.debug("Error?: stepHalfSize (%s) > signalSize (%s)",
                      stepHalfSize, signalSize)
        return np.zeros(signalSize, dtype=np.float_)

    result = np.zeros(signalSize, dtype=np.float_)
    if weight is not None:
        # Init weight sums
        highWeightSum = weight[:stepHalfSize].sum()
        # highSquareSum = np.exp2(weight[:stepHalfSize]).sum()
        highNonNormed = (weight[:stepHalfSize] * signal[:stepHalfSize]).sum()
        # Circular padding
        lowWeightSum = highWeightSum
        # lowSquareSum = highSquareSum
        lowNonNormed = -highNonNormed

    # ENH: vectorize this loop (it's the performance hotspot)
    for k in range(1, signalSize):
        highEnd = k + stepHalfSize - 1
        if highEnd >= signalSize:
            highEnd = signalSize - 1 - (highEnd - signalSize)
        lowEnd = k - stepHalfSize - 1
        if lowEnd < 0:
            lowEnd = -lowEnd - 1

        if weight is None:
            result[k] = result[k-1] + signal[highEnd] + signal[lowEnd] - 2*signal[k-1]
        else:
            lowNonNormed += signal[lowEnd] * weight[lowEnd] - signal[k-1] * weight[k-1]
            highNonNormed += signal[highEnd] * weight[highEnd] - signal[k-1] * weight[k-1]
            lowWeightSum += weight[k-1] - weight[lowEnd]
            highWeightSum += weight[highEnd] - weight[k-1]
            # lowSquareSum += weight[k-1] * weight[k-1] - weight[lowEnd] * weight[lowEnd]
            # highSquareSum += weight[highEnd] * weight[highEnd] - weight[k-1] * weight[k-1]
            result[k] = math.sqrt(stepHalfSize / 2) * (lowNonNormed / lowWeightSum +
                                                       highNonNormed / highWeightSum)

    if weight is None:
        stepNorm = math.sqrt(2. * stepHalfSize)
        result[1:signalSize] /= stepNorm

    return result


def FindLocalPeaks(signal, #const double * signal,
                   # peakLoc, #int * peakLoc
                  ):
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
        sig_prev, sig_curr, sig_next = signal[k-1:k+2]
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


def UnifyLevels(baseLevel, #const int * baseLevel,
                addonLevel, #const int * addonLevel,
                windowSize, #int windowSize,
                # joinedLevel, #int * joinedLevel);
               ):
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

    # Merge all addon items outside a window around each base item
    # ENH: do something clever with searchsorted & masks?
    joinedLevel = []
    addon_idx = 0
    for base_elem in baseLevel:
        while addon_idx < len(addonLevel):
            addon_elem = addonLevel[addon_idx]
            if addon_elem < base_elem - windowSize:
                # Addon is well before this base item -- use it
                joinedLevel.append(addon_elem)
                addon_idx += 1
            elif base_elem - windowSize <= addon_elem <= base_elem + windowSize:
                # Addon is too close to this base item -- skip it
                addon_idx += 1
            else:
                assert base_elem + windowSize < addon_elem
                # Addon is well beyond this base item -- keep for the next round
                break
        joinedLevel.append(base_elem)

    # Append the remaining addon items beyond the last base item's window
    last_pos = (baseLevel[-1] + windowSize if len(baseLevel) else -1)
    while addon_idx < len(addonLevel) and addonLevel[addon_idx] <= last_pos:
        addon_idx += 1
    if addon_idx < len(addonLevel):
        joinedLevel.extend(addonLevel[addon_idx:])

    return np.array(sorted(joinedLevel), dtype=np.int_)


def PulseConv(signal, #const double * signal,
              pulseSize, #int pulseSize,
             ):
    """Convolve a pulse function with a signal, applying circular padding to the
    signal.

    Used for non-stationary variance compensation.

    Parameters
    ----------
    signal: const array of floats
    pulseSize: int

    Returns
    -------
    array of floats

    Source: HaarSeg.c
    """
    signalSize = len(signal)
    if pulseSize > signalSize:
        # ENH: handle this endcase
        raise ValueError("pulseSize (%s) > signalSize (%s)"
                         % (pulseSize, signalSize))
    pulseHeight = 1. / pulseSize

    # Circular padding init
    result = np.zeros(signalSize, dtype=np.float_)
    for k in range((pulseSize + 1) // 2):
        result[0] += signal[k]
    for k in range(pulseSize // 2):
        result[0] += signal[k]
    result[0] *= pulseHeight

    n = 1
    for k in range(pulseSize // 2,
                   signalSize + (pulseSize // 2) - 1):
        tail = k - pulseSize
        if tail < 0:
            tail = -tail - 1
        head = k
        if head >= signalSize:
            head = signalSize - 1 - (head - signalSize)
        result[n] = result[n-1] + ((signal[head] - signal[tail]) * pulseHeight)
        n += 1

    return result


# XXX Apply afterward to the segmentation result? (not currently used)
def AdjustBreaks(signal, #const double * signal,
                 peakLoc, #const int * peakLoc,
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
        n1 = (npl_k if k == 0
                else npl_k - newPeakLoc[k-1])
        n2 = (len(signal) if k+1 == len(newPeakLoc)
              else newPeakLoc[k+1])- npl_k

        # Find the best offset for current breakpoint, trying only 1 sample
        # offset
        bestScore = float("Inf")  # Smaller is better
        bestOffset = 0
        for p in (-1, 0, 1):
            # Pointless to try to remove single-sample segments
            if (n1 == 1 and p == -1) or (n2 == 1 and p == 1):
                continue

            signal_n1_to_p = signal[npl_k - n1:npl_k + p]
            s1 = signal_n1_to_p.sum() / (n1 + p)
            ss1 = ((signal_n1_to_p - s1)**2).sum()

            signal_p_to_n2 = signal[npl_k + p:npl_k + n2]
            s2 = signal_p_to_n2.sum() / (n2 - p)
            ss2 = ((signal_p_to_n2 - s2)**2).sum()

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


if __name__ == '__main__':
    real_data = np.concatenate((np.zeros(800), np.ones(200),
                                np.zeros(800), .8*np.ones(200), np.zeros(800)))
    # np.random.seed(0x5EED)
    noisy_data = real_data + np.random.standard_normal(len(real_data)) * .2

    # # Run using default parameters
    seg_table = haarSeg(noisy_data, .005)

    logging.info("%s", seg_table)

    from matplotlib import pyplot
    indices = np.arange(len(noisy_data))
    pyplot.scatter(indices, noisy_data, alpha=0.2, color='gray')
    x, y = table2coords(seg_table)
    pyplot.plot(x, y, color='r', marker='x', lw=2)
    pyplot.show()

    # # The complete segmented signal
    # lines(seg.data$Segmented, col="red", lwd=3)
