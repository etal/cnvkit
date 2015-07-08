"""Call copy number variants from segmented log2 ratios."""

import numpy as np


def round_clonal(log2_ratio, purity, ploidy):
    """Infer integer copy number from log2 ratio.
    """
    return int(round(convert_diploid(log2_ratio)))


def round_diploid(log2_ratio):
    """Assume purity=1, ploidy=2."""
    return int(round(convert_diploid(log2_ratio)))


def round_thresholds(log2_ratio, thresholds=(-1.1, -0.3, 0.2, 0.7)):
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
        LOSS(1) <  -0.3
        GAIN(3) >=  +0.2
        AMP(4)  >=  +0.7

    For germline samples, better accuracy could be achieved with:

        LOSS(1) <  -0.4
        GAIN(3) >=  +0.3

    """
    i = 0
    for i, thresh in enumerate(thresholds):
        if log2_ratio < thresh:
            return i
    else:
        return int(np.ceil(convert_diploid(log2_ratio)))
        # ENH: opt. ploidy, purity
        # return int(np.ceil(convert_clonal(log2_ratio, p, p)))


def convert_clonal(log2_ratio, purity, ploidy):
    """
    """
    # TODO take math from export...


def convert_diploid(log2_ratio):
    """Assume purity=1, ploidy=2."""
    return 2 ** (log2_ratio + 1)


# Test: convert_clonal(x, 1, 2) == convert_diploid(x)

# _____________________________________________________________________________
# Rescaling etc.

# TODO refactor
def rescale_copy_ratios(cnarr, purity=None, ploidy=2, round_to_integer=False,
                        is_sample_female=None, is_reference_male=True):
    """Rescale segment copy ratio values given a known tumor purity."""
    if purity and not 0.0 < purity <= 1.0:
        raise ValueError("Purity must be between 0 and 1.")

    chr_x = core.guess_chr_x(cnarr)
    chr_y = ('chrY' if chr_x.startswith('chr') else 'Y')
    if is_sample_female is None:
        is_sample_female = core.guess_xx(cnarr, is_reference_male, chr_x,
                                         verbose=False)
    absolutes = cna_absolutes(cnarr, ploidy, purity, is_reference_male,
                              is_sample_female)
    if round_to_integer:
        absolutes = np.round(absolutes)
    # Avoid a logarithm domain error
    absolutes = np.maximum(absolutes, 0.0001)
    abslog = np.log2(absolutes / float(ploidy))
    newcnarr = cnarr.copy()
    newcnarr["coverage"] = abslog
    # Adjust sex chromosomes to be relative to the reference
    if is_reference_male:
        newcnarr['coverage'][newcnarr.chromosome == chr_x] += 1.0
    newcnarr['coverage'][newcnarr.chromosome == chr_y] += 1.0
    return newcnarr


def cna_absolutes(cnarr, ploidy, purity, is_reference_male, is_sample_female):
    """Calculate absolute copy number values from segment or bin log2 ratios."""
    absolutes = np.zeros(len(cnarr), dtype=np.float_)
    for i, row in enumerate(cnarr):
        ref_copies, expect_copies = _reference_expect_copies(
            row["chromosome"], ploidy, is_sample_female, is_reference_male)
        absolutes[i] = _log2_ratio_to_absolute(
            row["coverage"], ref_copies, expect_copies, purity)
    return absolutes


def _reference_expect_copies(chrom, ploidy, is_sample_female, is_reference_male):
    """Determine the number copies of a chromosome expected and in reference.

    For sex chromosomes, these values may not be the same ploidy as the
    autosomes.

    Return a pair: number of copies in the reference and expected in the sample.
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


def _log2_ratio_to_absolute(log2_ratio, ref_copies, expect_copies, purity=None):
    """Transform a log2 ratio value to absolute linear scale.

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

    If purity adjustment is skipped (p=1; e.g. if THetA was run beforehand):

        n = r*2^v
    """
    if purity and purity < 1.0:
        ncopies = (ref_copies * 2**log2_ratio - expect_copies * (1 - purity)
                  ) / purity
    else:
        ncopies = ref_copies * 2 ** log2_ratio
    return ncopies

