Additional scripts
==================

refFlat2bed.py
    Generate a BED file of the genes or exons in the reference genome given in
    UCSC refFlat.txt format.
    This script can be used in case the original BED file of targeted intervals
    is unavailable.  Subsequent steps of the pipeline will remove probes that
    did not receive sufficient coverage, including those exons or genes that
    were not targeted by the sequencing library.  However, better results are
    expected from CNVkit if the true targeted intervals can be provided.

genome2access.py:
    Calculate the sequence-accessible coordinates in chromosomes from the given
    reference genome, treating long spans of 'N' characters as the inaccessible
    regions.

