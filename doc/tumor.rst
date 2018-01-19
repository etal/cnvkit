Tumor analysis
==============

CNVkit has been used most extensively on solid tumor samples sequenced with a
target panel or whole-exome sequencing protocol. Several options and approaches
are available to support this use case:

- If you have unpaired tumor samples, or no normal samples sequenced on the same
  platform, see the :ref:`reference` command for strategies.

- Use ``--drop-low-coverage`` to ignore bins with log2 normalized coverage
  values below -15.  Virtually all tumor samples, even cancer cell lines, are
  not completely homogeneous. Even in regions of homozygous deletion in the
  largest tumor-cell clonal population, some sequencing reads will be obtained
  from contaminating normal cells without the deletion. Therefore, extremely low
  log2 copy ratio values do not indicate homozygous deletions but failed
  sequencing or mapping in all cells regardless of copy number status at that
  site, which are not informative for copy number. This option in the
  :ref:`batch` command applies to segmentation; the option is also available in
  the :ref:`segment`, :ref:`metrics`, :ref:`segmetrics`, :ref:`genemetrics` and
  :doc:`heterogeneity` commands.

    - Why -15? The null log2 value substituted for bins with zero coverage is
      -20 (about 1 millionth the average bin's coverage), and the maximum
      positive shift that can be introduced by normalizing to the reference is 5
      (for bins with 1/32 the average coverage; bins below this are masked out
      by the reference). In a .cnr file, any bins with log2 value below -15 are
      probably based on dummy values corresponding to zero-coverage (perhaps
      unmappable) bins, and not real observations.

- The :ref:`batch` command does not directly output integer copy number calls
  (see :doc:`heterogeneity`). Instead, use the ``--ploidy`` and ``--purity``
  options in :ref:`call` to calculate copy number for each sample individually
  using known or estimated tumor-cell fractions. Also consider using ``--center
  median`` in highly aneuploid samples to shift the log2 value of true neutral
  regions closer to zero, as it may be slightly off initially.

- If SNV calls are available in VCF format, use the ``-v``/``--vcf`` option in
  the :ref:`call` and :ref:`scatter` commands to calculate or plot b-allele
  frequencies alongside each segment's total copy number or log2 ratio. These
  values reveal allelic imbalance and loss of heterozygosity (LOH), supporting
  and extending the inferred CNVs.
