Chromosomal sex
===============

CNVkit attempts to handle chromosomal sex correctly throughout the analysis
pipelines. Several commands automatically infer a given sample's chromosomal
sex from the relative copy number of the autosomes and chromosomes X and Y; the
status log messages will indicate when this is happening. In most cases the
inference can be skipped or overridden by using the ``-x``/``--sample-sex``
option.

The :ref:`sex` command runs and report's CNVkit's inference for one or more
given samples, and can be used on .cnn, .cnr or .cns files at any stage of
processing.

Reference sex
-------------

See :ref:`reference`

If you want copy number calls to be relative to a male reference with a single X
chromosome but dipoid autosomes, use the ``-y`` option everywhere.
Otherwise, X and all autosomes will be considered normally diploid. Chromosome Y
will be considered haploid in either case.


Chromosomal sex in calling absolute copy number
-----------------------------------------------

See :ref:`call`

Plots and sex chromosomes
-------------------------

:ref:`diagram` adjusts the sex chromosomes for sample and reference sex so
that gains and losses on chromosomes X and Y are highlighted and labeled
appropriately.

:ref:`scatter` and :ref:`heatmap` do not adjust the sex chromosomes for sample
or reference sex.

FAQ
---

Why does chromosome X/Y show a gain/loss?
`````````````````````````````````````````

The copy number status of sex chromosomes may be surprising if you are unsure
about the sex of the samples and reference used:

- Female samples normalized to a male reference will show a doubling of
  chromosome X (log2 value about +1.0) and complete loss of chromosome Y (log2
  value below -3.0, usually far below).
- Male samples normalized to a female reference will show a single-copy loss of
  chromosome X (log2 value about -1.0). The chromosome Y value should be near
  0.0, but the log2 values may be noisier and less reliable than on autosomes.

In the output of the :ref:`diagram`, :ref:`call`, and :ref:`export` commands,
the X or Y copy number may be wrong if the sex of the reference
(``-y``/``--male-reference``) or sample (``-x``) was not specified correctly. If
sample sex was not specified on the command line, check the command's logged
status messages to see if the sample's sex was guessed incorrectly.

After you've verified the above, the CNV might be real.

CNVkit is not detecting my sample's sex correctly. What can I do?
`````````````````````````````````````````````````````````````````

In lower-quality samples, particularly tumor samples analyzed without a robust
reference (see :doc:`tumor`), there may be many bins with no coverage which bias
the segment means. Try repeating the :ref:`segment` command with the
``--drop-low-coverage`` option if you did not do so originally.

See also: https://www.biostars.org/p/210080/
