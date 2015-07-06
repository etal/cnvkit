Calling gains and losses
========================

The relationship between the observed copy ratio and the true underlying copy
number depends on tumor cell fraction (purity), genome ploidy (which may be
heterogeneous in a tissue sample), and the size of the subclonal population
containing the CNA. Because of these ambiguities, CNVkit only reports the
estimated log2 copy ratio, and does not currently attempt a formal statistical
test for estimating integer copy number.

In a diploid genome, a single-copy gain in a perfectly pure, homogeneous sample
has a copy ratio of 3/2. In log2 scale, this is log2(3/2) = 0.585, and a
single-copy loss is log2(1/2) = -1.0.

In the :ref:`diagram` plot, for the sake of providing a clean visualization of
confidently called CNAs, the default threshold to label genes is 0.5.  This
threshold will tend to display gene amplifications, fully clonal single-copy
gains in fairly pure samples, most single-copy losses, and complete deletions.

When using the :ref:`gainloss` command, I recommend choosing a threshold to suit
your needs depending on your knowledge of the sample's purity, heterogeneity,
and likely features of interest. As a starting point, try 0.1 or 0.2 if you are
going to do your own filtering downstream, or 0.3 if not.

Future releases of CNVkit will add another layer of statistical testing to make
variant calls from the log2 copy ratio data.

