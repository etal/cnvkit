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
command. In particular, the `--ci` option estimates the confidence interval of
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
