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
the sensitivity of the segmentation result.  This function includes several
optional extentions, supporting the use of weights (also known as quality of
measurments) and raw measurments.  We recommend using both extentions where
possible, as it greatly improves the segmentation result.  Raw red / green
measurments are used to detect low value probes, which are more sensitive to
noise.

Examples:

    real.data = c(rep.int(0,2000),rep.int(1,100),rep.int(0,2000));
    noisy.data = real.data + rnorm(length(real.data),sd = 0.2);
    plot(noisy.data)

    # Using default parameters
    seg.data = haarSeg(noisy.data);

    # Segments result table: segment start index | segment size | segment value
    print(seg.data$SegmentsTable)

    # The complete segmented signal
    lines(seg.data$Segmented, col="red", lwd=3)

"""
from __future__ import absolute_import, division, print_function

import math

import numpy as np
import pandas as pd
from scipy import stats

# from .. import params
from ..ngfrills import echo
from ..cnary import CopyNumArray as CNA


def segment_haar(cnarr):
    """Do segmentation for CNVkit.

    Calculate copy number segmentation by HaarSeg
    (http://haarseg.r-forge.r-project.org/)
    Input: log2 coverage data in Nexus 'basic' format
    Output: the CBS data table

    """
    chrom_tables = []
    # Segment each chromosome individually
    # ENH - skip large gaps (segment chrom. arms separately)
    for chrom, subprobes in cnarr.by_chromosome():
        # echo(chrom, ':')  # DBG
        segtable = haarSeg(subprobes['log2'])
        chromtable = pd.DataFrame({
            'chromosome': chrom,
            'start': np.asarray(subprobes['start']).take(segtable['start']),
            'end': np.asarray(subprobes['end']
                             ).take(segtable['start']+segtable['size']-1),
            'gene': '.',
            'log2': segtable['log2'],
            'probes': segtable['size'],
        })
        # echo(chromtable)  # DBG
        chrom_tables.append(chromtable)
    result = pd.concat(chrom_tables)
    echo("haar: Found", len(result), "segments")
    segarr = cnarr.as_dataframe(result)
    segarr.sort_columns()
    return segarr


# ---- from HaarSeg R code -- the API ----

def haarSeg(I,
            W=None,
            rawI=None,
            breaksFdrQ=0.005,  # orig. .001
            haarStartLevel=1,
            haarEndLevel=5):
    r"""Perform segmentation according to the HaarSeg algorithm.

    Arguments:

    I
        A 1D array of log-ratio values, sorted according to their genomic
        location.
    W
        Weight matrix, corresponding to quality of measurement, with values
        :math:`1/(\sigma^2)`. Must have the same size as I.
    rawI
        The minimum between the raw test-sample and control-sample coverages
        (before applying log ratio, but after any background reduction and/or
        normalization).
        Used for the non-stationary variance compensation.
        Must have the same size as I.
    breaksFdrQ
        The FDR q parameter. This value should lie between 0 and 0.5. The
        smaller this value is, the less sensitive the segmentation result will
        be.
        For example, we will detect fewer segmentation breaks when using Q =
        1e-4, compared to when using Q = 1e-3.
        Common used values are 1e-2, 1e-3, 1e-4.
    haarStartLevel
        The detail subband from which we start to detect peaks. The higher this
        value is, the less sensitive we are to short segments. The default is
        value is 1, corresponding to segments of 2 probes.
    haarEndLevel
        The detail subband until which we use to detect peaks. The higher this
        value is, the more sensitive we are to large trends in the data. This
        value DOES NOT indicate the largest possible segment that can be
        detected.  The default is value is 5, corresponding to step of 32 probes
        in each direction.


    Returns:

    A tuple of two elements:
        1. Segments result table, a list of tuples:
            (segment start index, segment size, segment value)
        2. Segmented signal array (same size as I)

    Source: haarSeg.R
    """
    def med_abs_diff(diff_vals):
        """Median absolute deviation, with deviations given."""
        if len(diff_vals) == 0:
            return 0.
        return np.median(np.abs(diff_vals)) / 0.6745

    diffI = HaarConv(I, None, 1)
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

        if rawI:
            pulseSize = 2 * stepHalfSize
            convMask = (PulseConv(varMask, pulseSize) >= .5)
            sigmaEst = (1 - convMask) * peakSigmaEst + convMask * noisySigmaEst
            convRes /= sigmaEst
            peakSigmaEst = 1.

        T = FDRThres(convRes[peakLoc], breaksFdrQ, peakSigmaEst)
        addonPeaks = HardThreshold(convRes, T, peakLoc)
        breakpoints = UnifyLevels(breakpoints, addonPeaks, 2 ** (level - 1))

    # echo("Found", len(breakpoints), "breakpoints:", breakpoints)

    segs = SegmentByPeaks(I, breakpoints, W)

    segSt = np.insert(breakpoints, 0, 0)
    segEd = np.append(breakpoints, len(I))
    segSize = segEd - segSt
    segValues = segs[segSt]

    segment_table = pd.DataFrame.from_items([('start', segSt),
                                             ('size', segSize),
                                             ('log2', segValues)])
    return segment_table


def FDRThres(x, q, stdev):
    """False discovery rate (FDR) threshold.

    Source: FDRThres.R
    """
    M = len(x)
    if M < 2:
        T = 0
    else:
        m = np.arange(1, M+1) / M
        sortedX = sorted(map(abs, x), reverse=True)

        # p = 2*(1 - pnorm(sortedX, sd = sdev));
        p = 2 * (1 - stats.norm.cdf(sortedX, stdev)) # XXX stdev

        # Get the largest index for which p <= m*q
        k = np.nonzero(p <= m * q)[0]
        if len(k):
            T = sortedX[k[-1]]
        else:
            # echo("No passing p-values: min p=%.4g, min m=%.4g, q=%s"
            #      % (p[0], m[0], q))
            # T = sortedX[1] + 1e-16;  # ~= 2^-52, like MATLAB "eps"
            T = sortedX[0] + 1e-16

    return T


def SegmentByPeaks(data, peaks, weights=None):
    """Average the values of the probes within each segment.

    `data` : the probe array values
    `peaks` : array of copy number breakpoints

    Source: SegmentByPeaks.R
    """
    segs = np.zeros_like(data)
    for seg_start, seg_end in zip(np.insert(peaks, 0, 0),
                                  np.append(peaks, len(data))):
        if weights is not None:
            # Weighted mean of individual probe values
            val = np.average(data[seg_start:seg_end],
                             weights=weights[seg_start:seg_end])
        else:
            # Unweighted mean of individual probe values
            val = np.mean(data[seg_start:seg_end])
        # echo("Segment value @%d-%d: %.4f" % (seg_start, seg_end, val))
        segs[seg_start:seg_end] = val
    return segs



# ---- from HaarSeg C code -- the core ----

# --- HaarSeg.h
def HaarConv(signal, #const double * signal,
             weight, #const double * weight,
             stepHalfSize, #int stepHalfSize,
            ):
    """Convolve haar wavelet function with a signal, applying circular padding
    to the signal.

    Supports weights when weight pointer is not None.

    Params:

        signal: const array of floats
        weight: const array of floats
        stepHalfSize: int

    Output:
        result: array of floats

    Source: HaarSeg.c
    """
    signalSize = len(signal)
    if stepHalfSize > signalSize:
        # XXX TODO handle this endcase
        # raise ValueError("stepHalfSize (%s) > signalSize (%s)"
        #                  % (stepHalfSize, signalSize))
        # echo("Error?: stepHalfSize (%s) > signalSize (%s)"
        #      % (stepHalfSize, signalSize))
        return np.zeros(signalSize, dtype=np.float_)

    result = np.zeros(signalSize, dtype=np.float_)
    if weight is not None:
        # Init weight sums
        highWeightSum = 0.
        highSquareSum = 0.
        highNonNormed = 0.
        for k in range(stepHalfSize):
            highWeightSum += weight[k]
            highSquareSum += weight[k]*weight[k]
            highNonNormed += weight[k]*signal[k]
        # Circular padding
        lowWeightSum = highWeightSum
        lowSquareSum = highSquareSum
        lowNonNormed = -highNonNormed

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
            lowSquareSum += weight[k-1] * weight[k-1] - weight[lowEnd] * weight[lowEnd]
            highSquareSum += weight[highEnd] * weight[highEnd] - weight[k-1] * weight[k-1]
            result[k] = math.sqrt(stepHalfSize / 2) * (
                        lowNonNormed / lowWeightSum + highNonNormed / highWeightSum)

    if weight is None:
        stepNorm = math.sqrt(2. * stepHalfSize)
        for k in range(1, signalSize):
            result[k] /= stepNorm

    return result


def FindLocalPeaks(signal, #const double * signal,
                   # peakLoc, #int * peakLoc
                  ):
    """Find local maxima on positive values, local minima on negative values.

    First and last index are never considered extramum.

    Parameters:

        signal: const array of floats

    Output:
        peakLoc: array of ints

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


def HardThreshold(signal, #const double * signal,
                  threshold, #double threshold,
                  peakLoc, #int * peakLoc);
                 ):
    """Drop any values of peakLoc under the given threshold.

    I.e. keep only the peak values where the signal amplitude is large enough.

    Parameters:

        signal: const array of floats
        threshold: scalar float
        peakLoc: modifiable array of ints

    Source: HaarSeg.c
    """
    peakLocOut = np.extract(np.abs(signal.take(peakLoc)) >= threshold, peakLoc)
    # echo("Peaks passing threshold:", peakLocOut,
    #      "\nat:", signal.take(peakLocOut))
    return peakLocOut


def UnifyLevels(baseLevel, #const int * baseLevel,
                addonLevel, #const int * addonLevel,
                windowSize, #int windowSize,
                # joinedLevel, #int * joinedLevel);
               ):
    """Unify several decomposition levels.

    Merge the two lists of breakpoints, but drop addonLevel values that are too
    close to baseLevel values.

    Parameters:

        baseLevel: const array of ints
        addonLevel: const array of ints
        windowSize: int

    Output:

        joinedLevel: array of ints

    Source: HaarSeg.c
    """
    joinedLevel = []

    # Going over all base
    addon_iter = iter(addonLevel)
    for base_elem in baseLevel:
        for addon_elem in addon_iter:
            if base_elem - windowSize <= addon_elem <= base_elem + windowSize:
                continue
            joinedLevel.append(addon_elem)
            if addon_elem > base_elem + windowSize:
                break

        joinedLevel.append(base_elem)

    # Insert remaining indexes in addon to joined
    joinedLevel.extend(addon_iter)
    return np.array(sorted(joinedLevel), dtype=np.int_)


# For the R version only, not Matlab (???)
def PulseConv(signal, #const double * signal,
              pulseSize, #int pulseSize,
             ):
    """Convolve a pulse function with a signal, applying circular padding to the
    signal.

    Parameters:

        signal: const array of floats
        pulseSize: int

    Output:
        result: modifiable array of floats

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

    Parameters:

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


# R array functions in Python
def match(x, table, nomatch=None):
    """Returns a vector of the positions of (first) matches of `x` in `table`.

    ## The intersection of two sets can be defined via match():
    ## Simple version:
    ## intersect <- function(x, y) y[match(x, y, nomatch = 0)]
    intersect # the R function in base, slightly more careful

    >>> match(range(5), range(2, 11))
    array([None, None, 0, 1, 2], dtype=object)
    >>> match(range(3,5), range(2, 11))
    array([1, 2])
    """
    return np.asarray([(list(table).index(val) if val in table else nomatch)
                       for val in x])


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
    seg_table = haarSeg(noisy_data)

    echo(seg_table)

    from matplotlib import pyplot
    indices = np.arange(len(noisy_data))
    pyplot.scatter(indices, noisy_data, alpha=0.2, color='gray')
    x, y = table2coords(seg_table)
    pyplot.plot(x, y, color='r', marker='x', lw=2)
    pyplot.show()

    # # The complete segmented signal
    # lines(seg.data$Segmented, col="red", lwd=3)
