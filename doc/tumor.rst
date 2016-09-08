Tumor analysis
==============

Solid tumor samples: Use ``--drop-low-coverage`` in the :ref:`batch` and
:ref:`segment` commands. Virtually all tumor samples, even cancer cell lines,
are not completely homogeneous. Even in regions of homozygous deletion in the
largest tumor-cell clonal population, some sequencing reads will be obtained
from contaminating normal cells without the deletion.
Therefore, extremely low log2 copy ratio values (below -15) do not indicate
homozygous deletions but failed sequencing or mapping in all cells regardless
of copy number status at that site, which are not informative for copy number.
This option in the :ref:`batch` command applies to segmentation; the option is
also available in the :ref:`segment`, :ref:`metrics`, :ref:`segmetrics`,
:ref:`gainloss` and :doc:`heterogeneity` commands.

If you have unpaired tumor samples, or no normal samples sequenced on the
same platform, see the :ref:`reference` command for strategies.
