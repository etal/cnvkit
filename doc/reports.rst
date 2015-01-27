Text and tabular reports
========================

breaks
------

List the targeted genes in which a segmentation breakpoint occurs.

::

    cnvkit.py breaks Sample.cnr Sample.cns

This helps to identify genes in which (a) an unbalanced fusion or other
structural rearrangement breakpoint occured, or (b) CNV calling is
simply difficult due to an inconsistent copy number signal.


gainloss
--------

Identify targeted genes with copy number gain or loss above or below a
threshold.

::

    cnvkit.py gainloss Sample.cnr
    cnvkit.py gainloss Sample.cnr -s Sample.cns -t 0.4 -y -m 5

If segments are given, the log2 ratio value reported for each gene will be the
value of the segment covering the gene. Where more than one segment overlaps the
gene, i.e. if the gene contains a breakpoint, each segment's value will be
reported as a separate row for the same gene.

If segments are not given, the median of the log2 ratio values of the bins
within each gene will be reported as the gene's overall log2 ratio value. This
mode will not attempt to identify breakpoints within genes.

The threshold (``-t``) and minimum number of bins (``-m``) options are used to
control which genes are reported. For example, a threshold of .6 (the default)
will report single-copy gains and losses in a completely pure tumor sample (or
germline CNVs), but a lower threshold would be necessary to call somatic CNAs if
significant normal-cell contamination is present.
Some likely false positives can be eliminated by dropping CNVs that cover a
small number of bins (e.g. with ``-m 3``, genes where only 1 or 2 bins show copy
number change will not be reported), at the risk of missing some true positives.

Specify the reference gender (``-y`` if male) to ensure CNVs on the X and Y
chromosomes are reported correctly; otherwise, a large number of spurious gains
or losses on the sex chromosomes may be reported.

The output is a text table of tab-separated values, which is amenable to further
processing by scripts and standard Unix tools such as ``grep``, ``sort``,
``cut`` and ``awk``.


gender
------

Guess samples' gender from the relative coverage of chrX. The answer is printed
to standard output as "Female" or "Male".

::

    cnvkit.py gender *.cnn *.cnr *.cns

metrics
-------

Compute coverage deviations and other metrics for self-evaluation.

::

    cnvkit.py metrics Sample.cnr -s Sample.cns

Multiple samples can be processed together to produce a table::

    cnvkit.py metrics S1.cnr S2.cnr -s S1.cns S2.cns

To compare multiple probe sets to the same segment dataset (reusing the given
segments for each coverage dataset)::

    cnvkit.py metrics Sample.targetcoverage.cnn Sample.antitargetcoverage.cnn Sample.cnr -s Sample.cns

(Note that the order of arguments and options matters here, unlike the other
commands: everything after the -s flag is treated as a segment dataset.)

