Quick start
===========

If you would like to quickly try CNVkit without installing it, try one of the
open-source `workflows on Dockstore
<https://dockstore.org/search?entryType=workflows&search=cnvkit>` with your platform of
choice. Separately maintained apps are also available on DNAnexus and Seven Bridges.

To run CNVkit on your own machine, keep reading.


Install CNVkit
--------------

See installation instructions in the README here:

https://github.com/etal/cnvkit


Download the reference genome
-----------------------------

Go to the `UCSC Genome Bioinformatics <http://hgdownload.soe.ucsc.edu/downloads.html>`_
website and download:

1. Your species' reference genome sequence, in FASTA format [required]
2. Gene annotation database, via RefSeq or Ensembl, in BED or "RefFlat" format
   (e.g.  `refFlat.txt
   <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz>`_)
   [optional]

You probably already have the reference genome sequence. If your species' genome
is not available from UCSC, use whatever reference sequence you have. CNVkit
only requires that your reference genome sequence be in FASTA format.
Both the reference genome sequence and the annotation database must be single,
uncompressed files.

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

If your sequencing protocol is :doc:`WGS <nonhybrid>`, then you don't need a
"target" BED file, but you probably do still want refFlat.txt.


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
reference genome sequence, and (optionally)
:ref:`sequencing-accessible regions <access>` along with your BAM files to:

1. Create a pooled :ref:`reference` of per-bin copy number estimates from
   several normal samples; then
2. Use this reference in processing all tumor samples that were sequenced with
   the same platform and library prep.

All of these steps are automated with the :ref:`batch` command. Assuming normal
samples share the suffix "Normal.bam" and tumor samples "Tumor.bam", a complete
``batch`` command could be::

    cnvkit.py batch *Tumor.bam --normal *Normal.bam \
        --targets my_baits.bed --fasta hg19.fasta \
        --access data/access-5kb-mappable.hg19.bed \
        --output-reference my_reference.cnn --output-dir example/

See the built-in help message to see what these options do, and for additional
options::

    cnvkit.py batch -h

If you have no normal samples to use for the :ref:`reference`, you can create a
"flat" reference which assumes equal coverage in all bins by using the
``--normal/-n`` flag without specifying any additional BAM files::

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

.. note:: **Which BED file should I use?**

    - *target* vs. *bait* BED files: For hybrid capture, the targeted regions
      (or "primary targets") are the genomic regions your capture kit attempts
      to ensure are well covered, e.g.  exons of genes of interest. The baited
      regions (or "capture targets") are the genomic regions your kit actually
      captures, usually including about 50bp flanking either side of each
      target. Give CNVkit the bait/capture BED file, not the primary targets.
    - For :ref:`wgs`, use the ``batch --method wgs`` option and optionally give
      the genome's "access" file -- if not given, it will be calculated from the
      genome sequence FASTA file.
    - For :ref:`tas`, use the ``batch --method amplicon`` option and give the
      target BED file.

    See also: :doc:`nonhybrid`

Next steps
----------

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

Now, starting a project from scratch, you could follow any of these approaches:

- Run ``batch`` as above with all tumor/test and normal/control samples
  specified as they are, and hope for the best. (This should usually work fine.)
- *For the careful:* Run ``batch`` with just the normal samples specified as
  normal, yielding coverage .cnn files and a **pooled reference**. Inspect the
  coverages of all samples with the :ref:`metrics` command, eliminating any
  poor-quality samples and choosing a larger or smaller antitarget bin size if
  necessary. Build an updated pooled reference using :ref:`batch` or
  :ref:`coverage` and :ref:`reference` (see :doc:`pipeline`), coordinating your
  work in a `Makefile <https://en.wikipedia.org/wiki/Makefile>`_, Rakefile, or
  similar build tool.

    - See also: `Ten Simple Rules for Reproducible Computational Research
      <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003285>`_

- *For the power user:* Run ``batch`` with all samples specified as tumor
  samples, using ``-n`` by itself to build a **flat reference**, yielding
  coverages, copy ratios, segments and optionally plots for all samples, both
  tumor and normal. Inspect the "rough draft" outputs and determine an
  appropriate strategy to build and use a **pooled reference** to re-analyze the
  samples -- ideally coordinated with a build tool as above.
- Use a framework like `bcbio-nextgen <https://bcbio-nextgen.readthedocs.io/>`_
  to coordinate the complete sequencing data analysis pipeline.

See the command-line usage pages for additional
:doc:`visualization <plots>`,
:doc:`reporting <reports>` and
:doc:`import/export <importexport>` commands in CNVkit.

