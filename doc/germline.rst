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

- The ``--filter`` option in :ref:`call` can be used to reduce the number of
  false-positive segments returned. To use the ``ci`` (recommended) or ``sem``
  filters, first run each sample's segmented .cns file through :ref:`segmetrics`
  with the ``--ci`` option, which adds upper and lower confidence limits to the
  .cns output that ``call --filter ci`` can then use.

- The ``--drop-low-coverage`` option (see :doc:`tumor`) should not be used; it
  will typically remove germline deep deletions altogether, which is not
  desirable.

- For using CNVkit with whole-genome sequencing datasets, see :doc:`nonhybrid`.
