"""Call copy number variants from segmented log2 ratios."""
from __future__ import absolute_import, division, print_function

import logging

import numpy as np
import pandas as pd

from . import segfilters


def do_call(cnarr, variants=None, method="threshold", ploidy=2, purity=None,
            is_reference_male=False, is_sample_female=False, filters=None,
            thresholds=(-1.1, -0.25, 0.2, 0.7)):
    if method not in ("threshold", "clonal", "none"):
        raise ValueError("Argument `method` must be one of: clonal, threshold")

    outarr = cnarr.copy()
    if filters:
        # Apply any filters that use segmetrics but not cn fields
        for filt in ('ci', 'sem'):
            if filt in filters:
                logging.info("Applying filter '%s'", filt)
                outarr = getattr(segfilters, filt)(outarr)
                filters.remove(filt)

    if variants:
        outarr["baf"] = variants.baf_by_ranges(outarr)

    if purity and purity < 1.0:
        logging.info("Rescaling sample with purity %g, ploidy %d",
                     purity, ploidy)
        absolutes = absolute_clonal(outarr, ploidy, purity,
                                    is_reference_male, is_sample_female)
        # Recalculate sample log2 ratios after rescaling for purity
        outarr["log2"] = log2_ratios(outarr, absolutes, ploidy,
                                     is_reference_male)
        if variants:
            # Rescale b-allele frequencies for purity
            outarr["baf"] = rescale_baf(purity, outarr["baf"])
    elif method == "clonal":
        # Estimate absolute copy numbers from the original log2 values
        logging.info("Calling copy number with clonal ploidy %d", ploidy)
        absolutes = absolute_pure(outarr, ploidy, is_reference_male)

    if method == "threshold":
        # Apply cutoffs to either original or rescaled log2 values
        tokens = ["%g => %d" % (thr, i) for i, thr in enumerate(thresholds)]
        logging.info("Calling copy number with thresholds: %s",
                     ", ".join(tokens))
        absolutes = absolute_threshold(outarr, ploidy, thresholds,
                                       is_reference_male)

    if method != 'none':
        outarr['cn'] = absolutes.round().astype('int')
        if 'baf' in outarr:
            # Calculate major and minor allelic copy numbers (s.t. cn1 >= cn2)
            upper_baf = ((outarr['baf'] - .5).abs() + .5).fillna(1.0).values
            outarr['cn1'] = ((absolutes * upper_baf).round()
                             .clip(0, outarr['cn'])
                             .astype('int'))
            outarr['cn2'] = outarr['cn'] - outarr['cn1']
            is_null = (outarr['baf'].isnull() & (outarr['cn'] > 0))
            outarr[is_null, 'cn1'] = np.nan
            outarr[is_null, 'cn2'] = np.nan

    if filters:
        # Apply the remaining cn-based filters
        for filt in filters:
            logging.info("Applying filter '%s'", filt)
            outarr = getattr(segfilters, filt)(outarr)

    outarr.sort_columns()
    return outarr


def log2_ratios(cnarr, absolutes, ploidy, is_reference_male,
                min_abs_val=1e-3, round_to_int=False):
    """Convert absolute copy numbers to log2 ratios.

    Optionally round copy numbers to integers.

    Account for reference sex & ploidy of sex chromosomes.
    """
    # Round absolute copy numbers to integer values
    if round_to_int:
        absolutes = absolutes.round()
    # Avoid a logarithm domain error
    ratios = np.log2(np.maximum(absolutes / ploidy, min_abs_val))
    # Adjust sex chromosomes to be relative to the reference
    if is_reference_male:
        ratios[(cnarr.chromosome == cnarr._chr_x_label).values] += 1.0
    ratios[(cnarr.chromosome == cnarr._chr_y_label).values] += 1.0
    return ratios


def absolute_threshold(cnarr, ploidy, thresholds, is_reference_male):
    """Call integer copy number using hard thresholds for each level.

    Integer values are assigned for log2 ratio values less than each given
    threshold value in sequence, counting up from zero.
    Above the last threshold value, integer copy numbers are called assuming
    full purity, diploidy, and rounding up.

    Default thresholds follow this heuristic for calling CNAs in a tumor sample:
    For single-copy gains and losses, assume 50% tumor cell clonality (including
    normal cell contamination). Then::

        R> log2(2:6 / 4)
        -1.0  -0.4150375  0.0  0.3219281  0.5849625

    Allowing for random noise of +/- 0.1, the cutoffs are::

        DEL(0)  <  -1.1
        LOSS(1) <  -0.25
        GAIN(3) >=  +0.2
        AMP(4)  >=  +0.7

    For germline samples, better precision could be achieved with::

        LOSS(1) <  -0.4
        GAIN(3) >=  +0.3

    """
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    for idx, row in enumerate(cnarr):
        cnum = 0
        ref_copies = _reference_copies_pure(row.chromosome, ploidy,
                                            is_reference_male)
        for cnum, thresh in enumerate(thresholds):
            if row.log2 <= thresh:
                if ref_copies != ploidy:
                    cnum = int(cnum * ref_copies / ploidy)
                break
        else:
            cnum = int(np.ceil(_log2_ratio_to_absolute_pure(row.log2,
                                                            ref_copies)))
        absolutes[idx] = cnum
    return absolutes


def absolute_clonal(cnarr, ploidy, purity, is_reference_male, is_sample_female):
    """Calculate absolute copy number values from segment or bin log2 ratios."""
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    for i, row in enumerate(cnarr):
        # TODO by_chromosome to reduce number of calls to this
        ref_copies, expect_copies = _reference_expect_copies(
            row.chromosome, ploidy, is_sample_female, is_reference_male)
        absolutes[i] = _log2_ratio_to_absolute(
            row.log2, ref_copies, expect_copies, purity)
    return absolutes


def absolute_pure(cnarr, ploidy, is_reference_male):
    """Calculate absolute copy number values from segment or bin log2 ratios."""
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    for i, row in enumerate(cnarr):
        ref_copies = _reference_copies_pure(row.chromosome, ploidy,
                                            is_reference_male)
        absolutes[i] = _log2_ratio_to_absolute_pure(row.log2, ref_copies)
    return absolutes


def absolute_dataframe(cnarr, ploidy, purity, is_reference_male, is_sample_female):
    """Absolute, expected and reference copy number in a DataFrame."""
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    reference_copies = expect_copies = np.zeros(len(cnarr), dtype=np.int_)
    for i, row in enumerate(cnarr):
        ref_copies, exp_copies = _reference_expect_copies(
            row.chromosome, ploidy, is_sample_female, is_reference_male)
        reference_copies[i] = ref_copies
        expect_copies[i] = exp_copies
        absolutes[i] = _log2_ratio_to_absolute(
            row.log2, ref_copies, exp_copies, purity)
    return pd.DataFrame({'absolute': absolutes,
                         'reference': reference_copies,
                         'expect': expect_copies})


def absolute_expect(cnarr, ploidy, is_sample_female):
    """Absolute integer number of expected copies in each bin.

    I.e. the given ploidy for autosomes, and XY or XX sex chromosome counts
    according to the sample's specified chromosomal sex.
    """
    exp_copies = np.repeat(ploidy, len(cnarr))
    is_y = (cnarr.chromosome == cnarr._chr_y_label).values
    if is_sample_female:
        exp_copies[is_y] = 0
    else:
        is_x = (cnarr.chromosome == cnarr._chr_x_label).values
        exp_copies[is_x | is_y] = ploidy // 2
    return exp_copies


def absolute_reference(cnarr, ploidy, is_reference_male):
    """Absolute integer number of reference copies in each bin.

    I.e. the given ploidy for autosomes, 1 or 2 X according to the reference
    sex, and always 1 copy of Y.
    """
    ref_copies = np.repeat(ploidy, len(cnarr))
    is_x = (cnarr.chromosome == cnarr._chr_x_label).values
    is_y = (cnarr.chromosome == cnarr._chr_y_label).values
    if is_reference_male:
        ref_copies[is_x] = ploidy // 2
    ref_copies[is_y] = ploidy // 2
    return ref_copies


def _reference_expect_copies(chrom, ploidy, is_sample_female, is_reference_male):
    """Determine the number copies of a chromosome expected and in reference.

    For sex chromosomes, these values may not be the same ploidy as the
    autosomes. The "reference" number is the chromosome's ploidy in the
    CNVkit reference, while "expect" is the chromosome's neutral ploidy in the
    given sample, based on the specified sex of each. E.g., given a female
    sample and a male reference, on chromosome X the "reference" value is 1 but
    "expect" is 2.

    Returns
    -------
    tuple
        A pair of integers: number of copies in the reference, and expected in
        the sample.
    """
    chrom = chrom.lower()
    if chrom in ["chrx", "x"]:
        ref_copies = (ploidy // 2 if is_reference_male else ploidy)
        exp_copies = (ploidy if is_sample_female else ploidy // 2)
    elif chrom in ["chry", "y"]:
        ref_copies = ploidy // 2
        exp_copies = (0 if is_sample_female else ploidy // 2)
    else:
        ref_copies = exp_copies = ploidy
    return ref_copies, exp_copies


def _reference_copies_pure(chrom, ploidy, is_reference_male):
    """Determine the reference number of chromosome copies (pure sample).

    Returns
    -------
    int
        Number of copies in the reference.
    """
    chrom = chrom.lower()
    if chrom in ["chry", "y"] or (is_reference_male and chrom in ["chrx", "x"]):
        ref_copies = ploidy // 2
    else:
        ref_copies = ploidy
    return ref_copies


def _log2_ratio_to_absolute(log2_ratio, ref_copies, expect_copies, purity=None):
    """Transform a log2 ratio to absolute linear scale (for an impure sample).

    Does not round to an integer absolute value here.

    Math::

        log2_ratio = log2(ncopies / ploidy)
        2^log2_ratio = ncopies / ploidy
        ncopies = ploidy * 2^log2_ratio

    With rescaling for purity::

        let v = log2 ratio value, p = tumor purity,
            r = reference ploidy, x = expected ploidy,
            n = tumor ploidy ("ncopies" above);

        v = log_2(p*n/r + (1-p)*x/r)
        2^v = p*n/r + (1-p)*x/r
        n*p/r = 2^v - (1-p)*x/r
        n = (r*2^v - x*(1-p)) / p

    If purity adjustment is skipped (p=1; e.g. if germline or if scaling for
    heterogeneity was done beforehand)::

        n = r*2^v
    """
    if purity and purity < 1.0:
        ncopies = (ref_copies * 2**log2_ratio - expect_copies * (1 - purity)
                  ) / purity
    else:
        ncopies = _log2_ratio_to_absolute_pure(log2_ratio, ref_copies)
    return ncopies


def _log2_ratio_to_absolute_pure(log2_ratio, ref_copies):
    """Transform a log2 ratio to absolute linear scale (for a pure sample).

    Purity adjustment is skipped. This is appropriate if the sample is germline
    or if scaling for tumor heterogeneity was done beforehand.

    .. math :: n = r*2^v
    """
    ncopies = ref_copies * 2 ** log2_ratio
    return ncopies


def rescale_baf(purity, observed_baf, normal_baf=0.5):
    """Adjust B-allele frequencies for sample purity.

    Math::

        t_baf*purity + n_baf*(1-purity) = obs_baf
        obs_baf - n_baf * (1-purity) = t_baf * purity
        t_baf = (obs_baf - n_baf * (1-purity))/purity
    """
    # ENH: use normal_baf array if available
    tumor_baf = (observed_baf - normal_baf * (1-purity)) / purity
    # ENH: warn if tumor_baf < 0 -- purity estimate may be too low
    return tumor_baf
