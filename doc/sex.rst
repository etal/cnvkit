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

When building a reference, each control sample's sex is inferred separately from
its target and antitarget coverage. These can disagree: antitarget bins on
chromosome Y are often too sparse to distinguish a male's Y from a female's
absent Y, which can make a genuinely male sample look female in the antitargets.
When the two sources conflict, CNVkit no longer simply defers to the antitargets;
instead it trusts whichever source has the more decisive chromosome-X coverage
(targets, by default, for capture panels). The status log reports the conflict
and the chosen source, e.g. ``...looks like male in targets but female in
antitargets; preferring targets``. You can always bypass inference by passing the
sample sex explicitly with ``-x``/``--sample-sex``.

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

For surfacing chrX/chrY variants in exports, the sex-aware view lives in VCF
and in BED via ``export bed --show variant``; both filter against the
sex-aware expected copy number (e.g. 1 on chrX in a male). The BED default
``--show ploidy`` is sex-agnostic by design (so a male's neutral haploid
chrX appears as ``cn=1`` rows and a chrX duplication to 2 copies is omitted
as "default ploidy" -- the format FreeBayes ``--cnv-map`` expects). See
:ref:`export` for the full rationale.

PAR handling
------------

Some reference genomes (files) have hard-masked (i.e. replaced with Ns) the
pseudoautosomal regions on chromosome Y. This means that male samples have
doubled coverage on chromosome X within PAR1/2 (since chromosome Y cannot be
covered at all). This will bias the reference creation as well as the copy
number calling. For example, in a mixed run with male & female samples, the
male samples are biased towards a copy number gain in PAR1/2 since they do have
similar coverage compared to the female samples. Conversely, the female samples
seem to exhibit a copy number loss there (they should have doubled coverage but
haven't).

In order to avoid this, the option ``--diploid-parx-genome <genome-build>`` can
be used for all affected commands. The genome build should be "grch38" or
similar. This will cause cnvkit to treat the PAR1/2 on chromosome X as
effectively autosomal (i.e. the expected copy number is the same for females
and males, the sex prediction expects the same coverage on X in PAR1/2 for both
sexes, etc.)

When ``--male-reference`` (see above) is also used, the reference copy number
for PAR1/2 on X is still 2 (same as for the autosomes).

See also: https://en.wikipedia.org/wiki/Pseudoautosomal_region
See also: https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components

Plots and sex chromosomes
-------------------------

:ref:`diagram` adjusts the sex chromosomes for sample and reference sex so
that gains and losses on chromosomes X and Y are highlighted and labeled
appropriately.

:ref:`scatter` and :ref:`heatmap` do not adjust the sex chromosomes for sample
or reference sex.

How sex is inferred from coverage
---------------------------------

For each input sample, CNVkit computes the difference between the median
log2 ratio of chromosome X and the median of the autosomes (and the same
for chromosome Y) and asks, for each chromosome, how close that observed
difference sits to the male expected position versus the female expected
position. On chromosome X the two expected positions are determined by
``-y`` / ``--haploid-x-reference``: a male reference puts male at log2=0
and female at log2=+1, while the default (diploid X) reference puts male
at log2=-1 and female at log2=0. On chromosome Y the female expected
position is the "no reads" sentinel (``NULL_LOG2_COVERAGE = -20``) and the
male expected position is autosome-median, which makes the chrY check act
as a wide-margin presence/absence detector.

The two per-chromosome "maleness" ratios are combined with an AND-gate: a
sample is called male only when *both* chrX and chrY independently look
more male than female. Female is the safe default when there isn't
positive evidence on both axes -- this is why a sample whose chrY is
masked or stripped (e.g. WGS with ``--diploid-parx-genome``) is called
female on chrX evidence alone only when there's still a Y signal to
gate on. Each sample's status log line shows both ratios, e.g.::

    Relative log2 coverage of chrX=-1.05, chrY=-21.4
    (maleness chrX=21, chrY=0.09) --> assuming female

Two ratios above 1.0 are required for male; either below 1.0 keeps the
female default. The 1.0 boundary is the midpoint between the two expected
positions in log-space, so the decision is threshold-free up to that
geometry.

When no chromosome-Y data is present at all -- for example, a female
reference profile that excluded chrY from the panel, or a target list with
no chrY bins -- the AND-gate naturally collapses to a chrX-only check
(``is_male = chrx_male_lr > 1.0``). The same 1.0 midpoint boundary still
applies, and the strict inequality means a sample sitting exactly at the
chrX midpoint ties to female (the safe default).

Non-human and Roman-numeral genomes
-----------------------------------

CNVkit recognizes chromosomes named with arabic numerals (``1`` / ``chr1``)
and Roman numerals (``chrI`` / ``XVI``). This means yeast (*S. cerevisiae*)
and other organisms that use Roman-numeral chromosome names are supported
out of the box.

For genomes where no recognizable autosome naming convention is detected
(e.g. an in-progress assembly with scaffold-style names), the
``autosomes()`` selection falls back permissively to returning the entire
dataset and logs a warning. This is intentional: it is better to plot
everything (and trust the user) than to silently discard data.

Sex-chromosome handling is *only* activated when both an ``X``/``chrX`` and
a ``Y``/``chrY`` (or just an ``X``/``chrX`` in a clearly arabic-numeral
genome) are detected in the data. In yeast — where ``chrX`` is the 10th
chromosome by Roman numeral rather than a sex chromosome —
sex-chromosome inference is automatically disabled. The same applies to
ZW-sex species (birds, reptiles, butterflies) whose chrZ/chrW are not
currently recognized, and to any custom assembly that doesn't include an
X/Y pair. In all these cases:

* ``cnvkit.py sex`` reports the sample as ``Unknown`` with ``NA`` for
  both log-ratios, rather than silently presenting the safe female
  default as a positive call.
* Commands that consume the inference internally (e.g. :ref:`call`,
  :ref:`diagram`) fall back to the safe female default and the
  per-sample ``Relative log2 coverage ...`` status line is simply
  omitted; nothing alarming is logged.

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
