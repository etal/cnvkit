Calling copy number gains and losses
====================================

The relationship between the observed copy ratio and the true underlying copy
number depends on tumor cell fraction (purity), genome ploidy (which may be
heterogeneous in a tissue sample), and the size of the subclonal population
containing the CNA. Because of these ambiguities, CNVkit initially reports only
the estimated log2 copy ratio, and supports several approaches to deriving
absolute integer copy number or evaluating copy number gains and losses.

In a diploid genome, a single-copy gain in a perfectly pure, homogeneous sample
has a copy ratio of 3/2. In log2 scale, this is log2(3/2) = 0.585, and a
single-copy loss is log2(1/2) = -1.0.

In the :ref:`diagram` plot, for the sake of providing a clean visualization of
confidently called CNAs, the default threshold to label genes is 0.5.  This
threshold will tend to display gene amplifications, fully clonal single-copy
gains in fairly pure samples, most single-copy losses, and total (homozygous)
deletions.

When using the :ref:`genemetrics` command, choose a threshold to suit your needs
depending on your knowledge of the sample's purity, heterogeneity, and likely
features of interest. As a starting point, try 0.1 or 0.2 if you are going to
do your own filtering downstream, or 0.3 if not.

To evaluate the level of support for each segment, use the :ref:`segmetrics`
command. In particular, the ``--ci`` option estimates the confidence interval of
the segment mean: the difference between the upper and lower limits suggests the
reliability of the estimate, and the value of the upper or lower limit suggests
the minimum loss or gain value, respectively, supported by the bin-level log2
ratios in that segment.

The :ref:`call` command implements two simple methods to convert the log2
ratios in a segmented .cns file to absolute integer copy number values. Given
known or estimated tumor purity and ploidy values, this command can also adjust
for :doc:`tumor heterogeneity <heterogeneity>`.

A .cns file can be converted to BED or VCF format using the :ref:`export`
command. These output styles display the inferred absolute integer copy number
value of each segment. To first adjust for known purity and ploidy of the
sample, or apply log2 thresholds for integer calls, use :ref:`call` before
exporting.


Filtering and merging
---------------------

The ``call`` command provides two distinct mechanisms for cleaning up
segmentation results: **filtering** (``--filter``) and **merging**
(``--merge``).

**Filtering** (``--filter``) neutralizes statistically uncertain segments
and merges them with adjacent neutral segments. This produces a clean call set
suitable for VCF export:

- ``--filter ci``: Use per-segment confidence intervals (from :ref:`segmetrics`
  ``--ci``). Segments whose CI overlaps zero are set to neutral and merged with
  adjacent neutrals. Segments confidently above or below zero are preserved
  individually.
- ``--filter sem``: Like ``ci``, but constructs the confidence interval from the
  standard error of the mean (from :ref:`segmetrics` ``--sem``).
- ``--filter cn``: After copy number calling, merge adjacent segments that have
  the expected neutral copy number (e.g. 2 for autosomes, 1 for chrX in males).
  Non-neutral segments are preserved individually even if adjacent segments share
  the same CN.
- ``--filter ampdel``: Drop all segments that are not amplifications or
  deletions, keeping only extreme events.

**Merging** (``--merge``) addresses over-segmentation by combining adjacent
segments that likely represent the same event:

- ``--merge bic``: Use the Bayesian Information Criterion to decide whether
  adjacent segments should be merged based on their log2 values. This runs
  before copy number calling. See :ref:`segmetrics` for details on BIC.
- ``--merge cn``: Merge any adjacent segments with the same integer copy number.
  Unlike ``--filter cn`` which only merges neutral segments, this merges all
  same-CN segments (e.g. two adjacent cn=1 losses).

Filters and merges can be combined. The order of operations is:

1. Pre-CN filters (``ci``, ``sem``) — mark uncertain segments neutral
2. Pre-CN merges (``bic``) — principled merging on log2 values
3. Copy number calling (threshold/clonal)
4. Post-CN filters (``cn``, ``ampdel``) — merge neutrals, drop non-events
5. Post-CN merges (``cn``) — merge adjacent same-CN segments

Example usage::

    # Filter uncertain segments, then merge over-segmentation
    cnvkit.py call sample.cns --filter ci --merge bic

    # Aggressive cleanup: filter + merge by copy number
    cnvkit.py call sample.cns --filter cn --merge cn

    # Multiple filters
    cnvkit.py call sample.cns --filter ci --filter cn
