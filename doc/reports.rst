Text and tabular reports
========================

.. _breaks:

``breaks``
----------

List the targeted genes in which a segmentation breakpoint occurs.

::

    cnvkit.py breaks Sample.cnr Sample.cns

This helps to identify genes in which (a) an unbalanced fusion or other
structural rearrangement breakpoint occured, or (b) CNV calling is
simply difficult due to an inconsistent copy number signal.

The output is a text table of tab-separated values, which is amenable to further
processing by scripts and standard Unix tools such as ``grep``, ``sort``,
``cut`` and ``awk``.

Columns:

- *gene*, *chromosome* -- as in .cns (see :ref:`cnxformat`), the gene where
  the breakpoint occurs and the chromosome it lies on.
- *location* -- the ``end`` of the segment to the left of the breakpoint, and
  ``start`` of the segment to the right.
- *change* -- the difference in ``log2`` values between the adjacent segments.
- *probes_left*, *probes_right* -- the number of probes on each side of the
  breakpoint within the gene. (Not the same as the number of probes supporting
  each segment; just the portion within the gene.)

For example, to get a list of the names of genes that contain a possible copy
number breakpoint (e.g. unbalanced translocation)::

    cnvkit.py breaks Sample.cnr Sample.cns | cut -f1 | sort -u > gene-breaks.txt


.. _genemetrics:

``genemetrics``
---------------

Identify targeted genes with copy number gain or loss above or below a
threshold. (Formerly called ``gainloss``.)

::

    cnvkit.py genemetrics Sample.cnr
    cnvkit.py genemetrics Sample.cnr -s Sample.cns -t 0.4 -m 5 -y

The first four columns of output table show each targeted gene's name and its
genomic coordinates (based on the first and last bins with that label in the
original target BED file, and thus the .cnr file).

The remaining output columns have slightly different meaning depending on
whether or not segments were provided. Without segments (.cnr alone):

- *log2*: Weighted mean of log2 ratios of all the gene's bins, including any
  off-target intronic bins.
- *depth*: Weighted mean of un-normalized read depths across all this gene's
  bins.
- *weight*: Sum of this gene's bins' weights.
- *nbins*: The number of bins assigned to this gene.

With segments (``-s``):

- *log2*: The log2 ratio value of the segment covering the gene, i.e. weighted
  mean of all bins covered by the whole segment, not just this gene.
- *depth*, *weight*, *probes*: As above.
- *seg_weight*: The sum of the weights of the bins supporting the segment.
- *seg_probes*: The number of probes supporting the segment.

The ``-t``/``--threshold`` and ``-m``/``--min-probes`` options are used to
control which genes are reported:

- A threshold of .2 (the default) will report single-copy gains and losses in a
  completely pure tumor sample (or germline CNVs), but a lower threshold would
  be necessary to call somatic CNAs if significant normal-cell contamination is
  present.
- Some likely false positives can be eliminated by dropping CNVs that cover a
  small number of bins, at the risk of missing some true positives.
  With ``-m 3``, the default, genes where only 1 or 2 bins show copy number
  change will not be reported.
  This applies to the segment's bin count (*seg_probes*) if segments are
  provided with ``-s``, otherwise it's the gene's bin count (*nbins*).

Specify the reference X-chromosome ploidy (``-y`` if the same option was used
when constructing the reference) to ensure CNVs on the X chromosome are reported
correctly; otherwise, a large number of spurious gains or losses may be
reported.

.. note:: Where more than one segment overlaps the gene, i.e. if the gene
    contains a breakpoint, each segment's value will be reported as a separate
    row for the same gene. If a large-scale CNA covers multiple genes, each of
    those genes will be listed individually.

The output is a text table of tab-separated values, like that of :ref:`breaks`.
Continuing the Unix example, we can try ``genemetrics`` both with and without the
segment files, take the intersection of those as a list of "trusted" genes, and
visualize each of them with :ref:`scatter`::

    cnvkit.py genemetrics -y Sample.cnr -s Sample.cns | tail -n+2 | cut -f1 | sort > segment-genes.txt
    cnvkit.py genemetrics -y Sample.cnr | tail -n+2 | cut -f1 | sort > ratio-genes.txt
    comm -12 ratio-genes.txt segment-genes.txt > trusted-genes.txt
    for gene in `cat trusted-genes.txt`
    do
        cnvkit.py scatter -s Sample.cn{s,r} -g $gene -o Sample-$gene-scatter.pdf
    done

(The point is that it's possible.)


.. _sex:

``sex``
-------

Guess samples' chromosomal sex from the relative coverage of chromosomes X and Y.
A table of the sample name (derived from the filename), inferred sex (string
"Female" or "Male"), and log2 ratio value of chromosomes X and Y is printed.

::

    cnvkit.py sex *.cnn *.cnr *.cns
    cnvkit.py sex -y *.cnn *.cnr *.cns

If there is any confusion in specifying either the sex of the sample or the
construction of the reference copy number profile, you can check what happened
using the ``sex`` command.
If the reference and intermediate .cnn files are available (.targetcoverage.cnn
and .antitargetcoverage.cnn, which are created before most of CNVkit's
corrections), CNVkit can report the reference sex and the sample's relative coverage of
the X and Y chromosomes::

    cnvkit.py sex reference.cnn Sample.targetcoverage.cnn Sample.antitargetcoverage.cnn

The output looks like this, where columns are filename, inferred sex, and
ratio of chromosome X and Y log2 coverages relative to autosomes::

    cnv_reference.cnn	Female	-0.176	-1.061
    Sample.targetcoverage.cnn	Female	-0.0818	-12.471
    Sample.antitargetcoverage.cnn	Female	-0.265	-15.139

If the ``-y`` option was not specified when constructing the reference (e.g.
``cnvkit.py batch ...``), then you have a female reference, and in the final
plots you will see chrX with neutral copy number in female samples and around -1
log2 ratio in male samples.


.. _metrics:

``metrics``
-----------

Calculate the spread of bin-level copy ratios from the corresponding final
segments using several statistics.
These statistics help quantify how "noisy" a sample is and help to decide which
samples to exclude from an analysis, or to select normal samples for a reference
copy number profile.

For a single sample::

    cnvkit.py metrics Sample.cnr -s Sample.cns

(Note that the order of arguments and options matters here, unlike the other
commands: Everything after the ``-s`` flag is treated as a segment dataset.)

Multiple samples can be processed together to produce a table::

    cnvkit.py metrics S1.cnr S2.cnr -s S1.cns S2.cns
    cnvkit.py metrics *.cnr -s *.cns

Several bin-level log2 ratio estimates for a single sample, such as the
uncorrected on- and off-target coverages and the final bin-level log2 ratios,
can be compared to the same final segmentation (reusing the given segments for
each coverage dataset)::

    cnvkit.py metrics Sample.targetcoverage.cnn Sample.antitargetcoverage.cnn Sample.cnr -s Sample.cns


In each case, given the bin-level copy ratios (.cnr) and segments (.cns) for a
sample, the log2 ratio value of each segment is subtracted from each of the bins
it covers, and several estimators of `spread
<https://en.wikipedia.org/wiki/Statistical_dispersion>`_ are calculated from the
residual values.
The output table shows for each sample:

- Total number of segments (in the .cns file) -- a large number of segments can
  indicate that the sample has either many real CNAs, or noisy coverage and
  therefore many spurious segments.
- Uncorrected sample `standard deviation
  <https://en.wikipedia.org/wiki/Standard_deviation>`_ -- this measure is prone
  to being inflated by a few outliers, such as may occur in regions of poor
  coverage or if the targets used with CNVkit analysis did not exactly match the
  capture. (Also note that the log2 ratio data are not quite normally
  distributed.) However, if a sample's standard deviation is drastically higher
  than the other estimates shown by the ``metrics`` command, that helpfully
  indicates the sample has some outlier bins.
- `Median absolute deviation
  <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_ (MAD) -- very
  `robust <https://en.wikipedia.org/wiki/Robust_measures_of_scale>`_ against
  outliers, but less `statistically efficient
  <https://en.wikipedia.org/wiki/Efficiency_%28statistics%29>`_.
- `Interquartile range <https://en.wikipedia.org/wiki/Interquartile_range>`_
  (IQR) -- another robust measure that is easy to understand.
- Tukey's `biweight midvariance
  <http://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidv.htm>`_
  -- a robust and efficient measure of spread.

Note that many small segments will fit noisy data better, shrinking the
residuals used to calculate the other estimates of spread, even if many of the
segments are spurious. One possible heuristic for judging the overall noisiness
of each sample in a table is to multiply the number of segments by the biweight
midvariance -- the value will tend to be higher for unreliable samples.
Check questionable samples for poor coverage (using e.g. `bedtools
<http://bedtools.readthedocs.io/>`_, `chanjo <http://www.chanjo.co/>`_,
`IGV <http://www.broadinstitute.org/igv/>`_ or `Picard CalculateHsMetrics
<http://broadinstitute.github.io/picard/command-line-overview.html#CalculateHsMetrics>`_).

Finally, visualizing a sample with CNVkit's :ref:`scatter` command will often
make it apparent whether a sample or the copy ratios within a genomic region can
be trusted.


.. _segmetrics:

``segmetrics``
--------------

Calculate summary statistics of the residual bin-level log2 ratio estimates
from the segment means, similar to the existing :ref:`metrics` command, but for each
segment individually.

Results are output in the same format as the CNVkit segmentation file (.cns),
with the stat names and calculated values printed in additional columns.

::

    cnvkit.py segmetrics Sample.cnr -s Sample.cns --iqr
    cnvkit.py segmetrics -s Sample.cn{s,r} --ci --pi

Supported stats:

- Alternative estimators of segment mean, which ignore bin weights: ``--mean``,
  ``-median``, ``--mode``.

- As in :ref:`metrics`: standard deviation (``--std``), median absolute
  deviation (``--mad``), inter-quartile range (``--iqr``), Tukey's biweight
  midvariance (``--bivar``)

- Additionally: mean squared error (``--mse``), standard error of the mean
  (``-sem``).

- Confidence interval of the segment mean (``--ci``), estimated by bootstrap
  (100 resamplings) of the bin-level log2 ratio values within the segment. The
  upper and lower bounds are output as separate columns ``ci_lo`` and ``ci_hi``.

- Prediction interval (``--pi``), estimated by the range between the 2.5-97.5
  percentiles of the segment's bin-level log2 ratios. The upper and lower bounds
  are output as columns ``pi_lo`` and ``pi_hi``.

The ``--ci`` and ``--sem`` values obtained here can also be used in the
:ref:`call` command for filtering segments.
