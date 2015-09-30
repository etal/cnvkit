"""Call copy number variants from segmented log2 ratios."""
from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd


def round_log2_ratios(cnarr, absolutes, ploidy, is_reference_male,
                      min_abs_val=1e-3):
    """Convert absolute integer copy numbers to log2 ratios.

    Account for reference gender & ploidy of sex chromosomes.
    """
    newcnarr = cnarr.copy()
    chr_x = cnarr._chr_x_label
    chr_y = cnarr._chr_y_label

    # Round absolute copy numbers to integer values
    absolutes = np.round(absolutes)
    # Avoid a logarithm domain error
    absolutes = np.maximum(absolutes, min_abs_val)
    newcnarr['log2'] = np.log2(absolutes / float(ploidy))

    # Adjust sex chromosomes to be relative to the reference
    if is_reference_male:
        newcnarr[newcnarr.chromosome == chr_x, 'log2'] += 1.0
    newcnarr[newcnarr.chromosome == chr_y, 'log2'] += 1.0
    return newcnarr


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
        ref_copies = _reference_copies_pure(row['chromosome'], ploidy,
                                            is_reference_male)
        for cnum, thresh in enumerate(thresholds):
            if row['log2'] <= thresh:
                if ref_copies != ploidy:
                    cnum = int(cnum * ref_copies / ploidy)
                break
        else:
            cnum = int(np.ceil(_log2_ratio_to_absolute_pure(row['log2'],
                                                            ref_copies)))
        absolutes[idx] = cnum
    return absolutes


def absolute_clonal(cnarr, ploidy, purity, is_reference_male, is_sample_female):
    """Calculate absolute copy number values from segment or bin log2 ratios."""
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    for i, row in enumerate(cnarr):
        ref_copies, expect_copies = _reference_expect_copies(
            row['chromosome'], ploidy, is_sample_female, is_reference_male)
        absolutes[i] = _log2_ratio_to_absolute(
            row['log2'], ref_copies, expect_copies, purity)
    return absolutes


def absolute_pure(cnarr, ploidy, is_reference_male):
    """Calculate absolute copy number values from segment or bin log2 ratios."""
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    for i, row in enumerate(cnarr):
        ref_copies = _reference_copies_pure(row['chromosome'], ploidy,
                                            is_reference_male)
        absolutes[i] = _log2_ratio_to_absolute_pure(row['log2'], ref_copies)
    return absolutes


def absolute_dataframe(cnarr, ploidy, purity, is_reference_male, is_sample_female):
    """Absolute, expected and reference copy number in a DataFrame."""
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    reference_copies = expect_copies = np.zeros(len(cnarr), dtype=np.int_)
    for i, row in enumerate(cnarr):
        ref_copies, exp_copies = _reference_expect_copies(
            row['chromosome'], ploidy, is_sample_female, is_reference_male)
        reference_copies[i] = ref_copies
        expect_copies[i] = exp_copies
        absolutes[i] = _log2_ratio_to_absolute(
            row['log2'], ref_copies, exp_copies, purity)
    return pd.DataFrame({'absolute': absolutes,
                         'reference': reference_copies,
                         'expect': expect_copies})


def _reference_expect_copies(chrom, ploidy, is_sample_female, is_reference_male):
    """Determine the number copies of a chromosome expected and in reference.

    For sex chromosomes, these values may not be the same ploidy as the
    autosomes. The "reference" number is the chromosome's ploidy in the
    CNVkit reference, while "expect" is the chromosome's neutral ploidy in the
    given sample, based on the specified gender of each. E.g., given a female
    sample and a male reference, on chromosome X the "reference" value is 1 but
    "expect" is 2.

    Return a pair of integers: number of copies in the reference, and expected in
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

    Return the integer number of copies in the reference.
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

    Math:

        log2_ratio = log2(ncopies / ploidy)
        2^log2_ratio = ncopies / ploidy
        ncopies = ploidy * 2^log2_ratio

    With rescaling for purity:

        let v = log2 ratio value, p = tumor purity,
            r = reference ploidy, x = expected ploidy;
        v = log_2(p*n/r + (1-p)*x/r)
        2^v = p*n/r + (1-p)*x/r
        n*p/r = 2^v - (1-p)*x/r
        n = (r*2^v - x*(1-p)) / p

    If purity adjustment is skipped (p=1; e.g. if germline or if scaling for
    heterogeneity was done beforehand):

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

    Math::

        n = r*2^v
    """
    ncopies = ref_copies * 2 ** log2_ratio
    return ncopies
