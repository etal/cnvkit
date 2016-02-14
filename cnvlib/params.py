"""Hard-coded parameters for CNVkit. These should not change between runs."""
# Filters used in constructing the reference
MIN_REF_COVERAGE = -5.0
MAX_REF_SPREAD = 1.0
NULL_LOG2_COVERAGE = -20

# Sequencing read length (or average, if it varies)
READ_LEN = 100

# Fragment size for paired-end reads
INSERT_SIZE = 250

# Target/bin names that are not meaningful gene names
# (In some UCSF panels, "CGH" probes denote selected intergenic regions)
IGNORE_GENE_NAMES = ("-", ".", "CGH")
