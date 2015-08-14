"""Robust estimators of central tendency and scale.

For use in evaluating performance of copy number estimation.

See:
    https://en.wikipedia.org/wiki/Robust_measures_of_scale
    https://astropy.readthedocs.org/en/latest/_modules/astropy/stats/funcs.html

"""
from __future__ import division

import numpy as np


def probe_deviations_from_segments(probes, segments, skip_low=True):
    """Difference in CN estimate of each probe from its segment."""
    probes.sort()
    segments.sort()
    if skip_low:
        probes = probes.drop_low_coverage()
    deviations = []
    for segment, subprobes in probes.by_segment(segments):
        deviations.append(subprobes['log2'] - segment['log2'])
    return np.concatenate(deviations)


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


def segment_mean(cnarr):
    """Weighted average of bin log2 values, ignoring too-low-coverage bins."""
    cnarr = cnarr.drop_low_coverage()
    if len(cnarr) == 0:
        return None
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
    https://astropy.readthedocs.org/en/latest/_modules/astropy/stats/funcs.html
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

