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

PAR1_X_GRCH37_START = 60001
PAR1_X_GRCH37_END = 2699520
PAR2_X_GRCH37_START = 154931044
PAR2_X_GRCH37_END = 155260560

PAR1_X_GRCH38_START = 10001
PAR1_X_GRCH38_END = 2781479
PAR2_X_GRCH38_START = 155701383
PAR2_X_GRCH38_END = 156030895