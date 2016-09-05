Quick start
===========

If you would like to quickly try CNVkit without installing it, try our app on
`DNAnexus <https://platform.dnanexus.com/app/cnvkit_batch>`_.

To run CNVkit on your own machine, keep reading.


Install CNVkit
--------------

Download the source code from GitHub:

https://github.com/etal/cnvkit

And read the README file.


Download the reference genome
-----------------------------

Go to the `UCSC Genome Bioinformatics <http://hgdownload.soe.ucsc.edu/downloads.html>`_
website and download:

1. Your species' reference genome sequence, in FASTA format [required]
2. Gene annotation database, via RefSeq or Ensembl, in "flat" format (e.g.
   `refFlat.txt
   <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz>`_)
   [optional]

You probably already have the reference genome sequence. If your species' genome
is not available from UCSC, use whatever reference sequence you have. CNVkit
only requires that your reference genome sequence be in FASTA format.
Both the reference genome sequence and the annotation database must be single,
uncompressed files.

**Sequencing-accessible regions:**
If your reference genome is the UCSC human genome hg19, a BED file of the
sequencing-accessible regions is included in the CNVkit distribution as
``data/access-5kb-mappable.hg19.bed``.
If you're not using hg19, consider building the "access" file yourself from your
reference genome sequence (say, ``mm10.fasta``) using the :ref:`access`
command::

    cnvkit.py access mm10.fasta -s 10000 -o access-10kb.mm10.bed

We'll use this file in the next step to ensure off-target bins ("antitargets")
are allocated only in chromosomal regions that can be mapped.

**Gene annotations:**
The gene annotations file (refFlat.txt) is useful to apply gene names to your
baits BED file, if the BED file does not already have short, informative names
for each bait interval. This file can be used in the next step.

If your targets look like::

    chr1	1508981	1509154
    chr1	2407978	2408183
    chr1	2409866	2410095

Then you want refFlat.txt.

Otherwise, if they look like::

    chr1	1508981	1509154	SSU72
    chr1	2407978	2408183	PLCH2
    chr1	2409866	2410095	PLCH2

Then you don't need refFlat.txt.


Map sequencing reads to the reference genome
--------------------------------------------

If you haven't done so already, use a sequence mapping/alignment program such as
`BWA <http://bio-bwa.sourceforge.net/>`_ to map your sequencing reads to the
reference genome sequence.

You should now have one or BAM files corresponding to individual samples.


Build a reference from normal samples and infer tumor copy ratios
-----------------------------------------------------------------

Here we'll assume the BAM files are a collection of "tumor" and "normal"
samples, although germline disease samples can be used equally well in place of
tumor samples.

CNVkit uses the bait BED file (provided by the vendor of your capture kit),
reference genome sequence, and sequencing-accessible regions along with your BAM
files to:

1. Create a pooled reference of per-bin copy number estimates from several
   normal samples; then
2. Use this reference in processing all tumor samples that were sequenced with
   the same platform and library prep.

All of these steps are automated with the :ref:`batch` command. Assuming normal
samples share the suffix "Normal.bam" and tumor samples "Tumor.bam", a complete
command could be::

    cnvkit.py batch *Tumor.bam --normal *Normal.bam \
        --targets my_baits.bed --fasta hg19.fasta \
        --access data/access-5kb-mappable.hg19.bed \
        --output-reference my_reference.cnn --output-dir example/

See the built-in help message to see what these options do, and for additional
options::

    cnvkit.py batch -h

If you have no normal samples to use for the reference, you can create a "flat"
reference which assumes equal coverage in all bins by using the ``--normal/-n``
flag without specifying any additional BAM files::

    cnvkit.py batch *Tumor.bam -n -t my_baits.bed -f hg19.fasta \
        --access data/access-5kb-mappable.hg19.bed \
        --output-reference my_flat_reference.cnn -d example2/

In either case, you should run this command with the reference genome sequence
FASTA file to extract GC and RepeatMasker information for bias corrections,
which enables CNVkit to improve the copy ratio estimates even without a paired
normal sample.

If your targets are missing gene names, you can add them here with the
``--annotate`` argument::

    cnvkit.py batch *Tumor.bam -n *Normal.bam -t my_baits.bed -f hg19.fasta \
        --annotate refFlat.txt --access data/access-5kb-mappable.hg19.bed \
        --output-reference my_flat_reference.cnn -d example3/


Process more tumor samples
--------------------------

You can reuse the reference file you've previously constructed to extract copy
number information from additional tumor sample BAM files, without repeating the
steps above.
Assuming the new tumor samples share the suffix "Tumor.bam" (and let's also
spread the workload across all available CPUs with the ``-p`` option, and
generate some figures)::

    cnvkit.py batch *Tumor.bam -r my_reference.cnn -p 0 --scatter --diagram -d example4/

The coordinates of the target and antitarget bins, the gene names for the
targets, and the GC and RepeatMasker information for bias corrections are
automatically extracted from the reference .cnn file you've built.

See the command-line usage pages for additional
:doc:`visualization <plots>`,
:doc:`reporting <reports>` and
:doc:`import/export <importexport>` commands in CNVkit.


Common Questions
----------------

Which BED file should I use?
````````````````````````````

- *target* vs. *bait* BED files: For hybrid capture, the targeted regions (or
  "primary targets") are the genomic regions your capture kit attempts to ensure
  are well covered, e.g.  exons of genes of interest. The baited regions (or
  "capture targets") are the genomic regions your kit actually captures, usually
  including about 50bp flanking either side of each target. Give CNVkit the
  bait/capture BED file, not the primary targets.
- For WGS, use the `batch --method wgs` option and optionally give the genome's
  "access" file -- if not given, it will be calculated from the genome sequence
  FASTA file.
- For targeted amplicon sequencing (TAS), use the `batch --method amplicon`
  option and give the target BED file.

How does CNVkit handle sample gender? Why does chrX show a gain/loss?
`````````````````````````````````````````````````````````````````````

If you want copy number calls to be relative to a male reference with a single X
chromosome but dipoid autosomes, use the `-y` option everywhere.
Otherwise, X and all autosomes will be considered normally diploid. Chromosome Y
will be considered haploid in either case.

How should I deal with tumor samples?
`````````````````````````````````````

- Solid tumor samples: Use --drop-low-coverage in the `batch` and `segment`
    commands. Virtually all tumor samples, even cancer cell lines, are
    not completely homogeneous. Even in regions of homozygous deletion in
    the largest tumor-cell clonal population, some sequencing reads will
    be obtained from contaminating normal cells without the deletion.
    Therefore, extremely low log2 copy ratio values (below -15) do not
    indicate homozygous deletions but failed sequencing or mapping
    in all cells regardless of copy number status at that site, which are
    not informative for copy number.
    This option in the `batch` command applies to segmentation; the
    option is also available in the `segment`, `metrics`, `segmetrics`,
    `gainloss` and `heterogeneity` commands.

- If you have unpaired tumor samples, or no normal samples sequenced on the
  same platform, see the :ref:`reference` command for strategies.

..  How should I deal with germline clinical samples?
..  `````````````````````````````````````````````````
