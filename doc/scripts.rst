Additional scripts
==================

refFlat2bed.py
    Generate a BED file of the genes or exons in the reference genome given in
    UCSC refFlat.txt format. (Download the input file from `UCSC Genome
    Bioinformatics <http://hgdownload.soe.ucsc.edu/downloads.html>`_).

    This script can be used in case the original BED file of targeted intervals
    is unavailable. Subsequent steps of the pipeline will remove probes that
    did not receive sufficient coverage, including those exons or genes that
    were not targeted by the sequencing library.  However, CNVkit will give much
    better results if the true targeted intervals can be provided.

reference2targets.py
    Extract target and antitarget BED files from a CNVkit reference file.
    While the :ref:`batch` command does this step automatically when an existing
    reference is provided, you may find this standalone script useful to recover
    the target and antitarget BED files that match the reference if those BED
    files are missing or you're not sure which ones are correct.

    Alternatively, once you have a stable CNVkit reference for your platform,
    you can use this script to drop the "bad" bins from your target and
    antitarget BED files (and subsequently built references) to avoid
    unnecessarily calculating coverage in those bins during future runs.

