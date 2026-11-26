"""Hard-coded parameters for CNVkit, and the small transforms that define them.

The values should not change between runs.
"""

import numpy as np
from numpy.typing import ArrayLike

# Filter thresholds used in constructing the reference (log2 scale)
MIN_REF_COVERAGE = -5.0
MAX_REF_SPREAD = 1.0
# The log2 value recorded for a bin with zero read depth, which has no finite
# logarithm. It is a placeholder, not a measurement: `log2_depth` assigns it,
# and `CopyNumArray.low_coverage_mask` recognizes it downstream.
NULL_LOG2_COVERAGE = -20.0


# Choices for the log2 value that `coverage` records for a bin with zero depth
ZERO_COVERAGE_FLOORS = ("fixed", "half-read")


def log2_depth(depth: ArrayLike, floor: ArrayLike = NULL_LOG2_COVERAGE) -> np.ndarray:
    """Convert read depths to log2 scale, with a hard floor at zero depth.

    A depth of zero becomes `floor`, as does a negative value, which is not a
    depth. Every positive value is transformed exactly, and a NaN depth stays
    NaN instead of being reported as zero coverage. Compare `safe_log2`, the
    soft floor of the RNA path, which also shifts the small nonzero values near
    its floor.

    Parameters
    ----------
    depth : array-like or float
        Absolute-scale read depths, or coverage ratios as Picard reports them.
    floor : array-like or float, optional
        The log2 value for a zero depth: one value, or one per element of
        `depth`. Default is `NULL_LOG2_COVERAGE`.

    Returns
    -------
    np.ndarray
        The log2-scale values, in the shape of `depth`.
    """
    depth = np.asarray(depth, dtype=float)
    log2 = np.full(depth.shape, floor, dtype=float)
    # Written so that NaN, for which every comparison is False, is transformed
    np.log2(depth, out=log2, where=~(depth <= 0))
    return log2


def half_read_log2(read_len: float, bin_sizes: ArrayLike) -> np.ndarray:
    """Get the log2 depth that half of one read would give each bin.

    This is a floor for bins that no read overlaps, in the units of the data:
    it lies one log2 unit below the depth of a single read aligned across the
    bin, at any bin size and sequencing depth. A read covers at most the width
    of a bin, so the floor of a bin narrower than a read is -1.

    Parameters
    ----------
    read_len : float
        Typical read length, in bases.
    bin_sizes : array-like
        Width of each bin, in bases.

    Returns
    -------
    np.ndarray
        The floor of each bin; `NULL_LOG2_COVERAGE` for a bin of no width, and
        for every bin if `read_len` is not a positive number.
    """
    bin_sizes = np.asarray(bin_sizes, dtype=float)
    floor = np.full(bin_sizes.shape, NULL_LOG2_COVERAGE)
    if read_len > 0:
        ok = bin_sizes > 0
        floor[ok] = np.log2(0.5 * np.minimum(read_len, bin_sizes[ok]) / bin_sizes[ok])
    return floor


def safe_log2(values, min_log2):
    """Transform values to log2 scale, with a soft floor at `min_log2`.

    Unlike `log2_depth`, this adds a constant to every value before the
    transform, so the result is smooth through zero but is not the exact log2
    of any input. The constant is absolute, so the values should already be
    normalized to a common scale, as the RNA path does to a median of 1.

    Parameters
    ----------
    values : np.array
        Absolute-scale values to transform. Should be non-negative.
    min_log2 : float
        Assign input zeros this log2-scaled value instead of -inf. Rather than
        hard-clipping, input values near 0 (especially below 2^min_log2) will be
        squeezed a bit above `min_log2` in the log2-scale output.
    """
    absolute_shift = 2**min_log2
    return np.log2(values + absolute_shift)


# Thresholds used in GC-content masking of bad bins at 'fix' step
GC_MIN_FRACTION = 0.3
GC_MAX_FRACTION = 0.7

# Fragment size for paired-end reads
INSERT_SIZE = 250

# Target/bin names that are not meaningful gene names
# (In some UCSF panels, "CGH" probes denote selected intergenic regions)
IGNORE_GENE_NAMES = ("-", ".", "CGH")
ANTITARGET_NAME = "Antitarget"
ANTITARGET_ALIASES = (ANTITARGET_NAME, "Background")

# PAR coordinates and the canonical names of supported genome builds are
# owned by skgenome.genomebuild; re-exported here for back-compat with code
# that imports them from cnvlib.params. Aliases such as "hg38" are accepted by
# skgenome.genomebuild.get_genome_build but are not keys here. Inner coordinate
# values are returned as mutable lists (matching the historical type) so any
# caller that treated them as such continues to work.
from skgenome.genomebuild import REGISTERED_BUILDS as _REGISTERED_BUILDS

PSEUDO_AUTSOMAL_REGIONS = {
    name: {region: list(coords) for region, coords in b.par_regions.items()}
    for name, b in _REGISTERED_BUILDS.items()
}
SUPPORTED_GENOMES_FOR_PAR_HANDLING = PSEUDO_AUTSOMAL_REGIONS.keys()
