"""Call copy number variants from segmented log2 ratios."""

from __future__ import annotations
import logging

import numpy as np

from . import segfilters
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cnvlib.cnary import CopyNumArray
    from cnvlib.vary import VariantArray
    from numpy import ndarray
    from pandas.core.frame import DataFrame
    from pandas.core.series import Series


def do_call(
    cnarr: CopyNumArray,
    variants: VariantArray | None = None,
    method: str = "threshold",
    ploidy: float = 2,
    purity: float | None = None,
    is_haploid_x_reference: bool = False,
    is_sample_female: bool = False,
    diploid_parx_genome: str | None = None,
    filters: list[str] | None = None,
    thresholds: tuple[float, float, float, float] | ndarray = (
        -1.1,
        -0.25,
        0.2,
        0.7,
    ),
) -> CopyNumArray:
    """Assign absolute integer copy number to each segment.

    This is the main API function for calling absolute copy numbers from
    segmented log2 ratios. It supports multiple calling methods and can
    optionally incorporate variant allele frequencies and tumor purity
    information.

    Parameters
    ----------
    cnarr : CopyNumArray
        Segmented copy number data (.cns file), typically from the segment
        command. Should contain 'log2' column with copy ratio values.
    variants : VariantArray, optional
        Variant allele frequency data from VCF, used to call allele-specific
        copy numbers. If provided, 'baf' (B-allele frequency) will be added
        to the output, and 'cn1'/'cn2' (major/minor allele copy numbers) will
        be calculated.
    method : str, optional
        Calling method to use. Options:
        - 'threshold': Apply log2 ratio thresholds (default)
        - 'clonal': Assume pure/clonal sample, infer from ploidy
        - 'none': Skip copy number calling, only apply filters
        Default: 'threshold'
    ploidy : float, optional
        Expected baseline ploidy of the sample. Usually 2 for diploid.
        Default: 2
    purity : float, optional
        Estimated tumor purity (0.0 to 1.0). If provided and < 1.0, log2
        ratios will be rescaled to account for normal cell contamination.
        Default: None (assume pure sample)
    is_haploid_x_reference : bool, optional
        Whether the reference sample is male (haploid X chromosome).
        Affects X chromosome copy number interpretation.
        Default: False
    is_sample_female : bool, optional
        Whether the test sample is female. Used with purity rescaling.
        Default: False
    diploid_parx_genome : str, optional
        Reference genome name for pseudo-autosomal region handling
        (e.g., 'hg19', 'hg38', 'mm10'). Treats PAR regions as diploid
        even in male samples.
        Default: None
    filters : list of str, optional
        Segment filters to apply. Options include 'ci' (confidence interval),
        'sem' (standard error), 'ampdel' (amplification/deletion), etc.
        See segfilters module for full list.
        Default: None
    thresholds : tuple of 4 floats or ndarray, optional
        Log2 ratio thresholds for calling copy numbers when method='threshold'.
        Format: (del_threshold, loss_threshold, gain_threshold, amp_threshold)
        These map to copy numbers [0, 1, 2, 3, 4+] respectively.
        Default: (-1.1, -0.25, 0.2, 0.7)

    Returns
    -------
    CopyNumArray
        Copy of input array with added 'cn' column (absolute copy number).
        If variants provided, also includes:
        - 'baf': B-allele frequency per segment
        - 'cn1': Major allele copy number (>= cn2)
        - 'cn2': Minor allele copy number (<= cn1)

    Raises
    ------
    ValueError
        If method is not one of 'threshold', 'clonal', or 'none'.

    See Also
    --------
    absolute_clonal : Calculate absolute copy numbers for pure samples
    absolute_threshold : Apply log2 thresholds to call copy numbers
    rescale_baf : Rescale B-allele frequencies for tumor purity

    Examples
    --------
    Basic threshold calling:
    >>> calls = do_call(segments, method='threshold')

    With tumor purity and variants:
    >>> calls = do_call(segments, variants=vcf, purity=0.7, ploidy=2)

    With custom thresholds:
    >>> calls = do_call(segments, thresholds=(-1.5, -0.3, 0.3, 1.0))
    """
    if method not in ("threshold", "clonal", "none"):
        raise ValueError("Argument `method` must be one of: clonal, threshold")

    outarr = cnarr.copy()
    if filters:
        # Apply any filters that use segmetrics but not cn fields
        for filt in ("ci", "sem"):
            if filt in filters:
                logging.info("Applying filter '%s'", filt)
                outarr = getattr(segfilters, filt)(outarr)
                filters.remove(filt)

    if variants:
        outarr["baf"] = variants.baf_by_ranges(outarr)

    if purity and purity < 1.0:
        logging.info("Rescaling sample with purity %g, ploidy %g", purity, ploidy)
        absolutes = absolute_clonal(
            outarr,
            ploidy,
            purity,
            is_haploid_x_reference,
            diploid_parx_genome,
            is_sample_female,
        )
        # Recalculate sample log2 ratios after rescaling for purity
        outarr["log2"] = log2_ratios(
            outarr, absolutes, ploidy, is_haploid_x_reference, diploid_parx_genome
        )
        if variants:
            # Rescale b-allele frequencies for purity
            outarr["baf"] = rescale_baf(purity, outarr["baf"])
    elif method == "clonal":
        # Estimate absolute copy numbers from the original log2 values
        logging.info("Calling copy number with clonal ploidy %g", ploidy)
        absolutes = absolute_pure(outarr, ploidy, is_haploid_x_reference)

    if method == "threshold":
        # Apply cutoffs to either original or rescaled log2 values
        tokens = ["%g => %d" % (thr, i) for i, thr in enumerate(thresholds)]
        logging.info("Calling copy number with thresholds: %s", ", ".join(tokens))
        absolutes = absolute_threshold(
            outarr, ploidy, thresholds, is_haploid_x_reference
        )

    if method != "none":
        outarr["cn"] = absolutes.round().astype("int")
        if "baf" in outarr:
            # Calculate major and minor allelic copy numbers (s.t. cn1 >= cn2)
            upper_baf = ((outarr["baf"] - 0.5).abs() + 0.5).fillna(1.0).to_numpy()
            outarr["cn1"] = (
                (absolutes * upper_baf).round().clip(0, outarr["cn"]).astype("int")
            )
            outarr["cn2"] = outarr["cn"] - outarr["cn1"]
            is_null = outarr["baf"].isna() & (outarr["cn"] > 0)
            outarr[is_null, "cn1"] = np.nan
            outarr[is_null, "cn2"] = np.nan

    if filters:
        # Apply the remaining cn-based filters
        for filt in filters:
            if not outarr.data.index.is_unique:
                logging.warning("Resetting index")  # DBG
                outarr.data = outarr.data.reset_index(drop=True)
            logging.warning("Applying filter '%s'", filt)
            outarr = getattr(segfilters, filt)(outarr)

    outarr.sort_columns()
    return outarr


def log2_ratios(
    cnarr: CopyNumArray,
    absolutes: Series,
    ploidy: float,
    is_haploid_x_reference: bool,
    diploid_parx_genome: str | None,
    min_abs_val: float = 1e-3,
    round_to_int: bool = False,
) -> Series:
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
    if is_haploid_x_reference:
        ratios[(cnarr.chr_x_filter(diploid_parx_genome)).to_numpy()] += 1.0
    ratios[(cnarr.chr_y_filter(diploid_parx_genome)).to_numpy()] += 1.0
    return ratios


def absolute_threshold(
    cnarr: CopyNumArray,
    ploidy: float,
    thresholds: tuple[float, float, float, float] | ndarray,
    is_haploid_x_reference: bool,
) -> ndarray:
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
    absolutes = np.zeros(len(cnarr), dtype=np.float64)
    for idx, row in enumerate(cnarr):
        ref_copies = _reference_copies_pure(
            row.chromosome, ploidy, is_haploid_x_reference
        )
        if np.isnan(row.log2):
            # XXX fallback
            logging.warning(
                "log2=nan found; replacing with neutral copy number %s", ref_copies
            )
            absolutes[idx] = ref_copies
            continue
        cnum = 0
        for cnum, thresh in enumerate(thresholds):
            if row.log2 <= thresh:
                if ref_copies != ploidy:
                    cnum = int(cnum * ref_copies / ploidy)
                break
        else:
            cnum = int(np.ceil(_log2_ratio_to_absolute_pure(row.log2, ref_copies)))
        absolutes[idx] = cnum
    return absolutes


def absolute_clonal(
    cnarr: CopyNumArray,
    ploidy: float,
    purity: float,
    is_haploid_x_reference: bool,
    diploid_parx_genome: str | None,
    is_sample_female: bool,
) -> Series:
    """Calculate absolute copy number values from segment or bin log2 ratios."""
    df = absolute_dataframe(
        cnarr,
        ploidy,
        purity,
        is_haploid_x_reference,
        diploid_parx_genome,
        is_sample_female,
    )

    return df["absolute"]


def absolute_pure(
    cnarr: CopyNumArray, ploidy: float, is_haploid_x_reference: bool
) -> ndarray:
    """Calculate absolute copy number values from segment or bin log2 ratios."""
    absolutes = np.zeros(len(cnarr), dtype=np.float64)
    for i, row in enumerate(cnarr):
        ref_copies = _reference_copies_pure(
            row.chromosome, ploidy, is_haploid_x_reference
        )
        absolutes[i] = _log2_ratio_to_absolute_pure(row.log2, ref_copies)
    return absolutes


def absolute_dataframe(
    cnarr: CopyNumArray,
    ploidy: float,
    purity: float,
    is_haploid_x_reference: bool,
    diploid_parx_genome: str | None,
    is_sample_female: bool,
) -> DataFrame:
    """Absolute, expected and reference copy number in a DataFrame."""
    df = get_as_dframe_and_set_reference_and_expect_copies(
        cnarr, ploidy, is_haploid_x_reference, diploid_parx_genome, is_sample_female
    )
    df["absolute"] = df.apply(
        lambda row: _log2_ratio_to_absolute(
            row["log2"], row["reference"], row["expect"], purity
        ),
        axis=1,
    )
    return df[["absolute", "expect", "reference"]]


def absolute_expect(
    cnarr: CopyNumArray,
    ploidy: float,
    diploid_parx_genome: str | None,
    is_sample_female: bool,
) -> Series:
    """Absolute integer number of expected copies in each bin.

    I.e. the given ploidy for autosomes, and XY or XX sex chromosome counts
    according to the sample's specified chromosomal sex.
    """
    is_haploid_x_reference = (
        True  # the reference sex doesn't matter for the expect column calculation
    )
    df = get_as_dframe_and_set_reference_and_expect_copies(
        cnarr, ploidy, is_haploid_x_reference, diploid_parx_genome, is_sample_female
    )
    exp_copies = df["expect"]
    return exp_copies


def absolute_reference(
    cnarr: CopyNumArray,
    ploidy: int,
    diploid_parx_genome: str | None,
    is_haploid_x_reference: bool,
) -> Series:
    """Absolute integer number of reference copies in each bin.

    I.e. the given ploidy for autosomes, 1 or 2 X according to the reference
    sex, and always 1 copy of Y.
    """
    is_sample_female = (
        True  # the sample sex doesn't matter for the reference column calculation
    )
    df = get_as_dframe_and_set_reference_and_expect_copies(
        cnarr, ploidy, is_haploid_x_reference, diploid_parx_genome, is_sample_female
    )
    ref_copies = df["reference"]
    return ref_copies


def get_as_dframe_and_set_reference_and_expect_copies(
    cnarr: CopyNumArray,
    ploidy: float,
    is_haploid_x_reference: bool,
    diploid_parx_genome: str | None,
    is_sample_female: bool,
) -> DataFrame:
    """Determine the number copies of a chromosome expected and in reference.

    For sex chromosomes, these values may not be the same ploidy as the
    autosomes. The "reference" number is the chromosome's ploidy in the
    CNVkit reference, while "expect" is the chromosome's neutral ploidy in the
    given sample, based on the specified sex of each. E.g., given a female
    sample and a male reference, on chromosome X the "reference" value is 1 but
    "expect" is 2. Note that the "reference" value for chromosome Y is always 1
    (better: `ploidy / 2`, see implementation below) to avoid divide-by-zero
    problems. The default reference is thus XXY (i.e. Klinefelter syndrome).

    Returns
    -------
    tuple
        A pair of integers: number of copies in the reference, and expected in
        the sample.
    """
    df = cnarr.copy().data

    # Set all to default value (i.e. for autsosomes):
    df["reference"] = np.repeat(ploidy, len(df))
    df["expect"] = np.repeat(ploidy, len(df))

    df.loc[cnarr.chr_x_filter(diploid_parx_genome), "reference"] = (
        ploidy // 2 if is_haploid_x_reference else ploidy
    )
    df.loc[cnarr.chr_x_filter(diploid_parx_genome), "expect"] = (
        ploidy if is_sample_female else ploidy // 2
    )

    df.loc[cnarr.chr_y_filter(diploid_parx_genome), "reference"] = ploidy // 2
    df.loc[cnarr.chr_y_filter(diploid_parx_genome), "expect"] = (
        0 if is_sample_female else ploidy // 2
    )
    if diploid_parx_genome is not None:
        # PAR1/2 are not covered on Y at all.
        df.loc[cnarr.pary_filter(diploid_parx_genome), "reference"] = 0
        df.loc[cnarr.pary_filter(diploid_parx_genome), "expect"] = 0
    return df


def _reference_copies_pure(
    chrom: str, ploidy: float, is_haploid_x_reference: bool
) -> float:
    """Determine the reference number of chromosome copies (pure sample).

    Returns
    -------
    float
        Number of copies in the reference.
    """
    chrom = chrom.lower()
    if chrom in ["chry", "y"] or (is_haploid_x_reference and chrom in ["chrx", "x"]):
        ref_copies = ploidy / 2
    else:
        ref_copies = ploidy
    return ref_copies


def _log2_ratio_to_absolute(
    log2_ratio: float,
    ref_copies: int,
    expect_copies: int,
    purity: float | None = None,
) -> float:
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
        ncopies = (ref_copies * 2**log2_ratio - expect_copies * (1 - purity)) / purity
    else:
        ncopies = _log2_ratio_to_absolute_pure(log2_ratio, ref_copies)
    return ncopies


def _log2_ratio_to_absolute_pure(log2_ratio: float, ref_copies: float) -> float:
    """Transform a log2 ratio to absolute linear scale (for a pure sample).

    Purity adjustment is skipped. This is appropriate if the sample is germline
    or if scaling for tumor heterogeneity was done beforehand.

    .. math :: n = r*2^v
    """
    ncopies = ref_copies * 2**log2_ratio
    return ncopies


def rescale_baf(purity: float, observed_baf: Series, normal_baf: float = 0.5) -> Series:
    """Adjust B-allele frequencies for sample purity.

    Math::

        t_baf*purity + n_baf*(1-purity) = obs_baf
        obs_baf - n_baf * (1-purity) = t_baf * purity
        t_baf = (obs_baf - n_baf * (1-purity))/purity
    """
    # ENH: use normal_baf array if available
    tumor_baf = (observed_baf - normal_baf * (1 - purity)) / purity
    # ENH: warn if tumor_baf < 0 -- purity estimate may be too low
    return tumor_baf
