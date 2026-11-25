"""Hard-coded parameters for CNVkit. These should not change between runs."""

import numpy as np

# Filter thresholds used in constructing the reference (log2 scale)
MIN_REF_COVERAGE = -5.0
MAX_REF_SPREAD = 1.0
# The log2 value recorded for a bin with zero read depth, which has no finite
# logarithm. It is a placeholder, not a measurement: `log2_depth` assigns it,
# and `CopyNumArray.low_coverage_mask` recognizes it downstream.
NULL_LOG2_COVERAGE = -20.0


def log2_depth(depth):
    """Convert read depths to log2 scale, with a hard floor at zero depth.

    A depth of exactly zero becomes `NULL_LOG2_COVERAGE`; every other value is
    transformed exactly, so a NaN depth stays NaN instead of being reported as
    zero coverage. Compare `safe_log2`, the soft floor of the RNA path, which
    also shifts the small nonzero values near its floor.

    Parameters
    ----------
    depth : array-like or float
        Absolute-scale read depths, or coverage ratios as Picard reports them.

    Returns
    -------
    np.ndarray
        The log2-scale values, in the shape of `depth`.
    """
    depth = np.asarray(depth, dtype=float)
    log2 = np.full(depth.shape, NULL_LOG2_COVERAGE)
    np.log2(depth, out=log2, where=depth != 0)
    return log2


def safe_log2(values, min_log2):
    """Transform values to log2 scale, with a soft floor at `min_log2`.

    Unlike `log2_depth`, this adds a constant to every value before the
    transform, so the result is smooth through zero but is exact nowhere. The
    constant is absolute, so the values should already be normalized to a
    common scale, as the RNA path does to a median of 1.

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
