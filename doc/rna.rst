RNA expression
==============

The relationship between genes' copy number and mRNA expression varies across
the genome. For a subset of genes, mostly housekeeping genes, the mRNA
expression levels measured by transcriptome sequencing are mostly explained by
underlying the genic regions' genomic copy number. CNVkit can use this
information to estimate coarse-grained copy number from RNA sequencing of a
process-matched cohort of samples.

Samples are processed simultaneously as a cohort, and two additional input files
are needed to complete the processing pipeline:

- Gene info -- a table of per-transcript and per-gene metadata exported from
  Ensembl BioMart. A snapshot of this file for the human transcriptome is
  bundled with CNVkit under ``data/ensembl-gene-info.hg38.tsv``.
- CNV-expression correlation coefficients -- calculated per gene via the script
  ``cnv_expression_correlate.py`` included with CNVkit. A pre-calculated set of
  coefficients calculated on TCGA melanoma datasets, obtained from cBioPortal,
  is bundled with CNVkit under ``data/tcga-skcm.cnv-expr-corr.tsv``.

With the above files, the RNA analysis can be run with either of the following
commands:

import-rsem
-----------

The third-party RNA quantification program RSEM produces
``*_rsem.genes.results`` output files that can be used as input to the
``import_rsem.py`` script included with CNVkit.

.. CNVkit's ``import-rsem`` command.

Use like this::

    import_rsem.py *_rsem.genes.results \
        --gene-resource data/ensembl-gene-info.hg38.tsv \
        --correlations data/tcga-skcm.cnv-expr-corr.tsv \
        --output out-summary.tsv --cnr-dir out/

Each gene's read counts and average transcript length are taken from each
sample's RSEM file. Normalized, bias-corrected ``*.cnr`` files are written to
``--cnr-dir``, and a summary table with all samples' data is written to
``-output``.


import-gene-counts
------------------

.. The ``import-gene-counts`` command 

The ``import_gene_counts.py`` script
can extract similar information from a 2-column, tab-delimited file of gene
names and per-gene read counts.

Use like this::

    import_gene_counts.py *.txt \
        --gene-resource data/ensembl-gene-info.hg38.tsv \
        --correlations data/tcga-skcm.cnv-expr-corr.tsv \
        --output out-summary.tsv --cnr-dir out/

Each gene's read counts and average transcript length are taken from each
sample's RSEM file. Normalized, bias-corrected ``*.cnr`` files are written to
``--cnr-dir``, and a summary table with all samples' data is written to
``-output``.


Segmentation
------------

The output ``*.cnr`` files can be used with CNVkit's :ref:`segment`,
:ref:`scatter`, and other commands similarly to the ``*.cnr`` files generated
from DNA sequencing data.  Each gene is represented by a single bin, and bins
may be overlapping or even nested.

The ``cbs`` segmentation method performs reasonably well to produce a
coarse-grained segmentation of the file. The ``hmm`` segmentation method also
does well.

When using the :ref:`scatter` command to  plot these files, note that the bin
(gene-level) weights are particularly important, and are visually represented by
circle size. A smoothed trendline (``-t``/``--trend``) can be helpful to
supplement the coarse-grained segmentation.


Considerations
--------------

- Input samples should be process-matched and ideally of the same source tissue
  type.
- The cohort size should be at least 5, preferably at least 10, samples.
- The ``--correlations`` input is not required but is strongly recommended. The
  TCGA melanoma cohort correlations can be used for analysis of any tissue type,
  not just neoplastic melanocytes. However, best results will usually be
  achieved with a correlations table specific to the test cohort. The script
  ``cnv_expression_correlate.py`` generates this table from input tables of
  per-gene and per-sample copy number and expression levels, typically retrieved
  from cBioPortal for TCGA cancer-specific cohorts.

