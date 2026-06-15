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

- **Gene info** -- a table of per-transcript and per-gene metadata exported from
  Ensembl BioMart. A snapshot of this file for the human transcriptome is
  bundled with CNVkit under ``data/ensembl-gene-info.hg38.tsv``.
- **CNV-expression correlation coefficients** -- calculated per gene via the
  script ``cnv_expression_correlate.py`` included with CNVkit. A pre-calculated
  set of coefficients calculated on TCGA melanoma datasets, obtained from
  cBioPortal, is bundled with CNVkit under ``data/tcga-skcm.cnv-expr-corr.tsv``.

With the above files, the RNA analysis can be run with either of the following
commands:

.. _import-rna:

``import-rna``
--------------


Use like this::

    cnvkit.py import-rna [ *_rsem.genes.results | *.txt ] \
        --gene-resource data/ensembl-gene-info.hg38.tsv \
        --correlations data/tcga-skcm.cnv-expr-corr.tsv \
        --output out-summary.tsv --output-dir out/

Each gene's read counts and average transcript length are taken from the input
file for each sample. Normalized, bias-corrected ``*.cnr`` files are written to
``--output-dir``, and an optional summary table with all samples' data is
written to ``-output``.

Input file sources:

- **RSEM:** The third-party RNA quantification program RSEM produces
  ``*_rsem.genes.results`` output files that can be used as input to the
  ``import-rna`` command.
- **Gene counts:** Alternatively, the gene Ensembl IDs and per-gene read counts
  can be read from a simple 2-column, tab-delimited file.
  This format is used by TCGA level 2 RNA expression data.
  You can also create the equivalent on your own from the output of another RNA
  quantification tool like Salmon or Kallisto.

Gene filtering
~~~~~~~~~~~~~~

Before normalization, ``import-rna`` discards genes that are rarely expressed
across the cohort: a gene is kept only if it has a detectable transcript (read
count of at least 1) in a minimum fraction of the input samples. By default this
fraction is 0.5, i.e. a gene must be expressed in at least half of the samples,
which reproduces the behavior of earlier CNVkit versions.

The ``--min-sample-fraction`` option lowers this threshold for single-cell or
otherwise sparse cohorts, where most genes are legitimately expressed in fewer
than half of the cells and the default would discard informative genes::

    cnvkit.py import-rna *.txt --gene-resource data/ensembl-gene-info.hg38.tsv \
        --min-sample-fraction 0.2 --output-dir out/

Lowering the threshold retains more genes in the output ``.cnr`` files, which can
shift downstream segmentation; this is a deliberate trade-off for sparse data.
Values below roughly 0.2 admit genes whose expression is detected in only a small
minority of samples, so their log2 ratios are likely dominated by noise.


Pseudo-autosomal regions on chromosome X
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pseudo-autosomal regions (PAR1 and PAR2) on chromosome X are diploid in both
sexes, so genes there behave like autosomal genes rather than sex-linked ones.
The ``--diploid-parx-genome`` option names a reference genome (e.g. ``grch38``)
whose PAR coordinates are then treated as autosomal when ``import-rna`` re-centers
each sample's log2 ratios::

    cnvkit.py import-rna *.txt --gene-resource data/ensembl-gene-info.hg38.tsv \
        --diploid-parx-genome grch38 --output-dir out/

Including the PAR-X bins in the autosomal pool shifts the centering value and so
changes the log2 output for the affected chromosome-X genes. When the option is
omitted (the default), output is unchanged.


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
  type. The DNA source for each sample can be single cells or bulk tissue.
- The cohort size should be at least 5, preferably at least 10, samples.
- The ``--correlations`` input is not required but is strongly recommended. The
  TCGA melanoma cohort correlations can be used for analysis of any tissue type,
  not just neoplastic melanocytes. However, best results will usually be
  achieved with a correlations table specific to the test cohort. The script
  ``cnv_expression_correlate.py`` generates this table from input tables of
  per-gene and per-sample copy number and expression levels, typically retrieved
  from cBioPortal for TCGA cancer-specific cohorts.


Using ``--normal`` to anchor against control samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``-n``/``--normal`` option accepts one or more count files representing
process-matched control (non-neoplastic) samples. When provided, the
normalization pipeline first applies a cohort-wide median polish to remove
per-sample sequencing-depth differences and per-gene scale differences, then
divides each gene's row by the median of the normal samples' values. Tumor
sample ``.cnr`` files are then expressed as log2 ratios against that
normal-anchored baseline.

This recentering has several important consequences that affect cohort design:

- **A single normal yields a non-informative ``.cnr`` for that sample.** With
  one normal, the per-gene "normal median" is just that sample's own value, so
  the anchoring step is a self-divide and the normal's output is approximately
  zero at every gene by construction. The normal's ``.cnr`` cannot be used to
  detect mislabelling, contamination, or copy-number events affecting the
  normal sample itself. ``import-rna`` emits a ``logging.warning`` in this
  case.
- **Two normals are still fragile.** A median over two values reduces to the
  arithmetic mean and remains a high-variance per-gene estimator. The output
  for each normal then represents that sample's deviation from the other
  normal, not a deviation from a stable biological baseline. ``import-rna``
  also warns in this case.
- **Three or more normals are recommended.** With ``n_normal >= 3``, the
  per-gene median is a robust estimator and individual normal ``.cnr`` files
  retain interpretable residuals (deviation from peer normals) usable for QC.
- **Cohort composition affects intermediate values.** The cohort-wide polish
  step uses the median across all samples in the cohort, so running the same
  ``--normal X`` on two different tumor cohorts produces different intermediate
  reference values. The final tumor log2 ratios are largely insensitive to
  this for typical cohorts, but extreme cohort imbalance (for example, a small
  number of tumors with shared, large-scale alterations) can perturb calls.
- **Control samples must themselves be free of the alterations under study.**
  Adjacent-tissue normals in solid-tumor settings frequently carry
  tumor-derived expression patterns (stromal contamination, field effects).
  When the normal cohort shares the same alterations as the tumors, the
  normal-anchored divide cancels real tumor signal and reduces the dynamic
  range of the resulting calls. Pure tissue-matched controls from independent
  donors are preferred when available.

If a suitable normal cohort is not available, run ``import-rna`` without
``--normal``. The output ``.cnr`` files will then be expressed as log2 ratios
against the cohort median, which carries its own assumption (that most genes
in most samples are copy-neutral) that may also be violated in heavily altered
cohorts.
