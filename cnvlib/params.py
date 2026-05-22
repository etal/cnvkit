"""Hard-coded parameters for CNVkit. These should not change between runs."""

# Filter thresholds used in constructing the reference (log2 scale)
MIN_REF_COVERAGE = -5.0
MAX_REF_SPREAD = 1.0
NULL_LOG2_COVERAGE = -20.0

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

# PAR coordinates and the list of supported genome builds are owned by
# skgenome.genomebuild; re-exported here for back-compat with code that
# imports them from cnvlib.params. Inner coordinate values are returned
# as mutable lists (matching the historical type) so any caller that
# treated them as such continues to work.
from skgenome.genomebuild import REGISTERED_BUILDS as _REGISTERED_BUILDS  # noqa: E402

PSEUDO_AUTSOMAL_REGIONS = {
    name: {region: list(coords) for region, coords in b.par_regions.items()}
    for name, b in _REGISTERED_BUILDS.items()
}
SUPPORTED_GENOMES_FOR_PAR_HANDLING = PSEUDO_AUTSOMAL_REGIONS.keys()
