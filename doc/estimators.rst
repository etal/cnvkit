Summary statistics and zero-coverage bins
=========================================

CNVkit reduces a set of bin-level log2 values to a single number in several
places: the log2 ratio of a segment, the log2 ratio of a gene in
:ref:`genemetrics`, the expected coverage of a bin across the samples pooled
into a :ref:`reference`, and the center that each sample is shifted to. These
sites do not all use the same estimator. This page states which estimator each
site uses and the reasoning that decides it.


Why the estimator matters
-------------------------

The log2 of a bin's read depth is approximately normally distributed where the
bin holds many reads. Under that condition the mean of the bins in a region is
the most efficient estimate of the region's log2 ratio, and weighting each bin
by its reliability improves the estimate further.

The approximation fails as the number of reads per bin falls. Read counts are
discrete and bounded at zero, so the log2 values of bins that hold only a few
reads are skewed toward low values, and a bin with no reads has no finite log2
value at all. CNVkit records such a bin at a placeholder log2 value of -20 in
the output of the :ref:`coverage` command. The placeholder is not a
measurement: it marks a value as censored, somewhere below what a single read
would have shown.

A mean has no resistance to a value of this kind. One bin at the placeholder
among *n* equally weighted bins lowers their mean by about 20/*n* log2 units,
so a single empty bin among eight copy-neutral bins yields -2.5, which is well
below the -1 of a single-copy loss in a pure sample. A robust estimator, such
as the median or Tukey's biweight location, is unaffected by the placeholder as
long as fewer than half of the values carry it. The same property makes it
unaffected by a real deletion that covers fewer than half of the bins, and it
is less efficient than the mean where every value is a measurement.

Neither kind of estimator is correct everywhere, so the choice follows from
what the number is meant to represent.


Which estimator each site uses
------------------------------

**A value that measures one genomic interval is the weighted mean of the log2
values of the bins in that interval**, with each bin weighted by its ``weight``
column. Bins removed before segmentation, for example by
``--drop-low-coverage``, do not contribute. This applies to:

- The ``log2`` value of every segment written by the :ref:`segment` command,
  whichever segmentation method produced the breakpoints.
- The ``log2`` value of a segment that the :ref:`call` command merges from
  several others, which is the weighted mean of the merged segments and
  therefore of their bins.
- The gene-level ``log2`` value that :ref:`genemetrics` reports when no
  segments are given.
- A gene that the :ref:`diagram` command draws as a feature of its own, so that
  the plot agrees with the :ref:`genemetrics` table.

One statistic at all of these sites means that a threshold has the same meaning
at each of them: the value of a gene and the value of the segment covering it
are the same quantity, computed over different sets of bins. The confidence
interval that :ref:`segmetrics` and :ref:`genemetrics` report with ``--ci``,
and that ``call --filter ci`` compares with zero, is a bootstrap interval for
this weighted mean.

**A value that summarizes a heterogeneous set, which has no single true value,
is a robust estimator.** This applies to:

- The expected log2 coverage of a bin in a pooled :ref:`reference`, which is
  Tukey's biweight location of that bin across the normal samples, with the
  biweight midvariance as its spread. One atypical sample does not move the
  expected value.
- The center that :ref:`fix`, :ref:`reference` and ``call --center`` shift a
  sample to, which by default is the median of the per-chromosome medians of
  the autosomal bins. Aneuploid chromosomes and unevenly targeted panels do not
  move it.
- The chromosome X and Y log2 ratios used to infer sex, which are medians
  relative to the autosomal median (see :doc:`sex`).
- The variance of a sample's bins that :ref:`fix` uses to set bin weights,
  which is the biweight midvariance of the bins' residuals, calculated without
  the zero-coverage bins.
- A :ref:`diagram` feature that combines many bins or genes sharing one drawn
  point, which is the biweight location of those bins.

The ``--mean``, ``--median`` and ``--mode`` options of :ref:`segmetrics` and
:ref:`genemetrics` add unweighted alternatives beside the weighted mean; they
do not replace it.


Zero-coverage bins
------------------

A bin is treated as having no usable coverage if its read depth is zero or
missing, or if its log2 value is below -15. The second condition identifies the
placeholder in a .cnr file, where normalization to the reference has shifted it
from -20 by at most 5 (see :doc:`tumor`), and in files that lack a ``depth``
column.

Each site handles these bins as follows:

- :ref:`reference`: a bin without coverage in most of the pooled samples
  receives a pooled log2 value below -5, and is then excluded from every
  sample normalized to that reference. Where only a minority of the pooled
  samples lack coverage in a bin, the biweight location disregards them.
- :ref:`fix`: the final centering of each sample and the variance used for bin
  weights are calculated without these bins. The bins themselves remain in the
  .cnr file.
- :ref:`segment`: the bins are segmented with all others by default, which
  lowers the mean of each segment that contains them. The command logs a
  warning when at least 1% of the autosomal bins are affected. The
  ``--drop-low-coverage`` option excludes them.
- :ref:`genemetrics`, :ref:`segmetrics` and :ref:`metrics`: the bins are
  included by default. The ``--drop-low-coverage`` option excludes them from the
  log2 value and from every statistic reported with it.

The bins are kept by default because they can be the signal. In a germline
sample or a tumor sample of high purity, a homozygous deletion leaves bins with
no reads, and excluding those bins would report the deleted gene as copy-neutral
within a longer segment, or omit it from the :ref:`genemetrics` table. In a
bulk tumor sample, where normal cells contribute reads to every deleted region,
a bin with no reads more often reflects a region that could not be sequenced or
mapped, and :doc:`tumor` recommends ``--drop-low-coverage`` for that case.

Where bins lack reads because they are too small for the sequencing depth, as
in shallow whole-genome sequencing, neither choice is satisfactory: the empty
bins are sampling zeros scattered among bins that hold one or two reads.
Recomputing coverage with larger bins (see :ref:`autobin`) restores the
condition under which the weighted mean is appropriate, without discarding
data.

The size of the placeholder's effect on a mean follows from its distance from
the data, which is arbitrary. The ``--zero-coverage-floor half-read`` option of
the :ref:`coverage` command (see :ref:`zero_coverage_floor`) records an empty
bin one log2 unit below a bin that holds a single read instead, so that its
effect on a mean is in proportion to the depth of the data. Such bins are
still recognized by their depth of 0, and each site above handles them in the
same way. A pooled :ref:`reference` still excludes a bin that most of the
pooled samples left empty, because the pooled depth of that bin is 0.


Reading gene-level values
-------------------------

Beside each gene's ``probes`` count, the :ref:`genemetrics` table reports
``low_probes``, the number of the gene's bins that are treated as having no
usable coverage. Comparing the two shows how far the placeholder accounts for
the gene's ``log2`` value:

- Where ``low_probes`` is 0, every bin of the gene is a measurement.
- Where ``low_probes`` equals ``probes``, no bin of the gene has coverage. This
  is the pattern of a homozygous deletion of the whole gene, and also that of
  a gene that cannot be sequenced or mapped in this assay.
- Where ``low_probes`` lies between the two, ``log2`` combines measurements
  with placeholder values, and its magnitude reflects the fraction of empty
  bins more than a copy number. This is the pattern of a partial deletion, and
  also that of a few bins that cannot be sequenced or mapped. The bin-level
  values in the .cnr file distinguish the two, and repeating the command with
  ``--drop-low-coverage`` reports the weighted mean of the remaining bins.

The ``--median`` option adds a value that does not respond to a minority of
the gene's bins. A ``log2`` far below the ``median`` indicates that a minority
of the bins are much lower than the rest, whether they are empty and counted
in ``low_probes`` or hold only a few reads and are not.

The ``--ci`` interval is calculated from the same bins as ``log2``, so an
interval that is wide relative to those of other genes with a similar number of
bins indicates the same disagreement among bins.
