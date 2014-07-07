Text and tabular reports
========================

breaks
------

List the targeted genes in which a segmentation breakpoint occurs.

::

    cnvkit.py breaks Sample.cnr Sample.cns

gainloss
--------

Identify targeted genes with copy number gain or loss.

::

    cnvkit.py gainloss Sample.cnr
    cnvkit.py gainloss Sample.cnr Sample.cns -t 0.4

gender
------

Guess samples' gender from the relative coverage of chrX. Answer is printed to
standard output.

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

