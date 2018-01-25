Chromosomal sex
===============

CNVkit attempts to handle chromosomal sex correctly at every step.
Several commands automatically infer a given sample's chromosomal sex from the
relative copy number of the autosomes and chromosomes X and Y; the status log
messages will indicate when this is happening.
In most cases the inference can be skipped or overridden by using the
``-x``/``--sample-sex`` option.

The :ref:`sex` command runs and report's CNVkit's inference for one or more
given samples, and can be used on .cnn, .cnr or .cns files at any stage of
processing.

Reference sex-chromosome ploidy
-------------------------------

By default, copy number calls and log2 ratios will be relative to a diploid X
chromosome and haploid Y. This happens regardless of the control samples used in
the :ref:`reference` command; if any input samples are male (haploid X), their X
chromosome coverage depths are doubled in order to be equivalent to diploid X.

However, some users prefer calls relative to a haploid X chromosome -- this is
how array CGH data are usually presented, for example. In that context it's
often referred to as a "male reference". This convention can be enabled in
CNVkit by using the ``-y`` / ``--male-reference`` / ``--haploid-x-reference``
option.  Note that this does not require any of the control samples to be male;
female control samples' X coverage depth is automatically halved so that it
appears as haploid in the CNVkit pipeline. Chromosome Y is always treated as
haploid in either case.

Chromosomal sex in calling absolute copy number
-----------------------------------------------

If ``-y`` is used to construct a reference, then the same option should also be
used with the commands :ref:`call`, :ref:`export`, and :ref:`genemetrics` for
samples processed with that reference.

Note that the options ``-x``/``--sample-sex`` and ``-y`` / ``--male-reference``
/ ``--haploid-x-reference`` are different: If a female sample is run with a
haploid-X reference, segments on chromosome X with log2-ratio +1 will be treated
as copy-number-neutral, because that's the expected copy number, while an
X-chromosome segment with log2-ratio 0 will be treated as a hemizygous loss.


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
the X-chromosome copy number may be wrong if the X-ploidy of the
reference (``-y``) or sample (``-x``) was not specified correctly. If
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
