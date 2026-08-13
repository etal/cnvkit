Germline analysis
=================

CNVkit can be used with exome sequencing of constitutional (non-tumor) samples,
for example to detect germline copy number alterations associated with heritable
conditions. However, note that CNVkit is less accurate in detecting CNVs
smaller than 1 Mbp, typically only detecting variants that span multiple exons
or captured regions.  When used on exome or target panel datasets, CNVkit will
not detect the small CNVs that are more common in populations.

To use CNVkit to detect medium-to-large CNVs or unbalanced SVs in constitutional
samples:

- The :ref:`call` command can be used directly without specifying the
  ``--purity`` and ``--ploidy`` values, as the defaults will be correct for
  mammalian cells. (For non-diploid species, use the correct ``--ploidy``, of
  course.) The default ``--method threshold`` assigns integer copy number
  similarly to ``--method clonal``, but with smaller thresholds for calling
  single-copy changes. The default thresholds allow for mosaicism in CNVs, which
  have smaller log2 value than a single-copy CNV would indicate. (They're more
  common than often thought.)

- The ``--filter`` option in :ref:`call` reduces the number of false-positive
  segments returned; ``ci`` is the recommended filter, and the example below
  shows the two commands it takes.

- The ``--drop-low-coverage`` option (see :doc:`tumor`) should not be used; it
  will typically remove germline deep deletions altogether, which is not
  desirable.


A worked example
----------------

Analyze the case samples against a pooled reference built from a set of control
samples captured and sequenced the same way. Assuming the case samples are
named ``Case1.bam``, ``Case2.bam`` and so on, and the controls
``Control1.bam`` and so on::

    cnvkit.py batch Case*.bam --normal Control*.bam \
        --targets my_baits.bed --annotate data/refFlat_hg38.txt \
        --fasta hg38.fasta --access data/access-10kb.hg38.bed \
        --output-reference my_reference.cnn --output-dir results/

The control samples contribute to the reference only. Each case sample yields
its own set of files in ``results/``, among them the segmented log2 ratios
``Case1.cns`` and the integer copy numbers ``Case1.call.cns``, assigned with
the default calling method and thresholds; the complete set of outputs is
listed under :ref:`batch`.

The ``ci`` filter that ``batch`` applies while calling uses a permissive
confidence level, corresponding to ``--alpha 0.5`` in :ref:`segmetrics`. To
filter the segments more strictly, or to use any calling option that ``batch``
does not expose, run :ref:`segmetrics` and :ref:`call` on a case sample's
segmented ``.cns`` file. The default ``--alpha 0.05`` yields a wider
confidence interval, so more of the uncertain segments are set to neutral and
merged with their neighbors::

    cnvkit.py segmetrics results/Case1.cnr -s results/Case1.cns --ci \
        -o results/Case1.segmetrics.cns
    cnvkit.py call results/Case1.segmetrics.cns --filter ci \
        -o results/Case1.germline.cns

Repeat both commands for each case sample.


Whole-genome and single-sample data
-----------------------------------

A single sample analyzed against a flat reference is usually too noisy to
support calling the small copy number variants that constitutional analysis is
concerned with. A pooled reference is strongly preferable here, as it is for
tumor samples and for the same reasons -- see :ref:`reference` -- even where
the contributing samples are not matched controls. For single-sample
whole-genome data in particular, a dedicated structural variant caller such as
Manta or Parliament2 is better suited to the task. For the whole-genome
workflow itself, see :doc:`nonhybrid`.
