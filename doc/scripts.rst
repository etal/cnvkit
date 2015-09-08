Additional scripts
==================

genome2access.py:
    Calculate the sequence-accessible coordinates in chromosomes from the given
    reference genome, treating long spans of 'N' characters as the inaccessible
    regions and outputting the coordinates of the regions between them.
    An "access" file for the UCSC hg19 genome is included in the CNVkit source
    distribution under the ``data/`` directory.

    CNVkit will compute "antitarget" bins only within the accessible genomic
    regions specified in the "access" file produced by this script. If there are
    many small excluded/inaccessible regions in the genome, then small,
    less-reliable antitarget bins would be squeezed into the remaining
    accessible regions.  The ``-s`` option tells the script to ignore short
    regions that would otherwise be excluded as inaccessible, allowing larger
    antitarget bins to overlap them.

    Additional regions to exclude can also be given with the ``-x`` option. This
    option can be used more than once to exclude several BED files listing
    different sets of regions. For example, "excludable" regions of poor
    mappability have been precalculated by others and are available from the
    `UCSC FTP Server <ftp://hgdownload.soe.ucsc.edu/goldenPath/>`_
    (see `here for hg19
    <ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/>`_).


refFlat2bed.py
    Generate a BED file of the genes or exons in the reference genome given in
    UCSC refFlat.txt format.  (Download the input file from `UCSC Genome
    Bioinformatics <http://hgdownload.soe.ucsc.edu/downloads.html>`_).

    This script can be used in case the original BED file of targeted intervals
    is unavailable. Subsequent steps of the pipeline will remove probes that
    did not receive sufficient coverage, including those exons or genes that
    were not targeted by the sequencing library.  However, CNVkit will give much
    better results if the true targeted intervals can be provided.

