Whole-genome sequencing and targeted amplicon capture
=====================================================

CNVkit is designed for use on **hybrid capture** sequencing data, where
off-target reads are present and can be used improve copy number estimates.

If necessary, CNVkit can be used on **whole-genome sequencing** (WGS) datasets
by specifying the genome's sequencing-accessible regions as the "targets",
avoiding "antitargets", and using a gene annotation database to label genes in
the resulting BED file::

    cnvkit.py batch ... -t data/access-10000.hg19.bed -g data/access-10000.hg19.bed --split --annotate refFlat.txt

Or::

    cnvkit.py target data/access-10000.hg19.bed --split --annotate refFlat.txt -o Targets.bed
    cnvkit.py antitarget data/access-10000.hg19.bed -g data/access-10000.hg19.bed -o Background.bed

This produces a "target" binning of the entire sequencing-accessible area of the
genome, and empty "antitarget" files which CNVkit will handle safely from
version 0.3.4 onward.


Similarly, to use CNVkit on **targeted amplicon sequencing** data instead --
although this is not recommended -- you can exclude all off-target regions from
the analysis by passing the target BED file as the "access" file as well::

    cnvkit.py batch ... -t Targeted.bed -g Targeted.bed ...

Or::

    cnvkit.py antitarget Targeted.bed -g Targeted.bed -o Background.bed

However, this approach does not collect any copy number information between
targeted regions, so it should only be used if you have in fact prepared your
samples with a targeted amplicon sequencing protocol.  It also does not attempt
to normalize each amplicon at the gene level, though this may be addressed in a
future version of CNVkit.

