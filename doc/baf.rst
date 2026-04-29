Allele frequencies and copy number
==================================

What is BAF?
------------

In this context, the "B" allele is the non-reference allele observed in a germline
heterozygous SNP, i.e. in the normal/control sample. Since the tumor cells' DNA
originally derived from normal cells' DNA, most of these SNPs will also be
present in the tumor sample. But due to allele-specific copy number alterations,
loss of heterozygosity or allelic imbalance, the allelic frequency of these SNPs
may be different in the tumor, and that's evidence that one (or both) of the
germline copies was gained or lost during tumor evolution.

The shift in b-allele frequency is calculated relative to the expected
heterozygous frequency 0.5, and minor allele frequencies are "mirrored" above
and below 0.5 so that it does not matter which allele is considered the
reference -- the relative shift from 0.5 will be the same either way. (Multiple
alternate alleles are not considered here.)


How does it work?
-----------------

Estimation from SNP b-allele frequencies works by comparing the shift in
allele frequency of heterozygous, germline SNPs in the tumor sample from the
expected ~50% -- e.g. a 3-copy segment in a diploid genome would have log2
ratio of +0.58 and heterozygous SNPs would have an average BAF of 67% or 33% if
the tumor sample is fully clonal, and closer to log2 0.0 and BAF 0.5 if there
is normal-cell contamination and/or :doc:`tumor heterogeneity <heterogeneity>`.

Use SNP b-allele frequencies from a VCF in these commands:

- :ref:`call`
- :ref:`scatter`
- :ref:`export` ``nexus-ogt`` and ``theta``


.. _baf-vcf-prep:

Preparing a VCF
---------------

CNVkit's BAF analysis requires a VCF that contains heterozygous **germline**
SNVs in the analyzed sample. Somatic-only variant calls cannot be used --
their allele frequencies do not encode the allelic-imbalance signal that BAF
analysis depends on. The clearest setup is:

- A joint tumor/normal VCF produced by an SNV caller such as FreeBayes,
  VarDict, or MuTect, with both samples present.
- Germline SNVs are kept in the file (do not strip them out). For MuTect,
  this typically means keeping the ``REJECT`` records.
- Somatic variants are flagged with the ``SOMATIC`` tag in the INFO column.
  CNVkit will skip these by default.

If you do not have a matched normal sample, you can use 1000 Genomes common
SNP sites to extract likely germline SNVs from a tumor-only VCF and use just
those sites. Note that without a paired normal, sample-level somatic filtering
is weaker (see below).

If you already have somatic calls produced by GATK Mutect2, you can re-run with
``--genotype-germline-sites true --genotype-pon-sites true`` to retain the
germline SNP sites in the same VCF -- this gives CNVkit something to work with.

An `example VCF
<https://github.com/etal/cnvkit/blob/master/test/formats/na12878_na12882_mix.vcf?raw=true>`_
constructed from the 1000 Genomes samples NA12878 and NA12882 is included in
CNVkit's test suite.

See also: :ref:`vcfformat`.


Sample identification in VCFs
-----------------------------

When the VCF contains multiple samples, CNVkit needs to know which is the
tumor and which is the normal. There are two ways:

- **PEDIGREE header tag.** Add a PEDIGREE tag to the VCF header declaring the
  tumor sample(s) as ``Derived`` and the normal as ``Original``. CNVkit will
  detect this automatically.
- **Command-line options.** Pass ``-i``/``--sample-id`` to identify the tumor
  sample and ``-n``/``--normal-id`` to identify the matched normal. These
  options are accepted by the :ref:`call`, :ref:`scatter`, :ref:`segment`,
  :ref:`export` ``theta``, and :ref:`export` ``nexus-ogt`` commands.

If neither is provided, CNVkit silently uses the first sample in the VCF and
treats the input as unpaired -- so genotype-based somatic filtering cannot
run, and a tumor-with-somatic-only VCF will look superficially heterozygous
to CNVkit. This is a common cause of incorrect BAF output (see
:ref:`baf-troubleshooting` below).


.. _baf-rescaling:

BAF rescaling for purity
------------------------

When ``cnvkit call --purity P`` is run with a VCF, the per-segment observed
BAF is rescaled to estimate the BAF of pure tumor cells. The model assumes
the observed BAF is a mixture of tumor and normal contributions:

.. math::

    \text{obs\_baf} = p \cdot \text{tumor\_baf} + (1 - p) \cdot \text{normal\_baf}

where :math:`p` is the tumor purity and the normal contribution is at the
heterozygous baseline :math:`\text{normal\_baf} = 0.5`. Solving for the
tumor BAF:

.. math::

    \text{tumor\_baf} = \frac{\text{obs\_baf} - 0.5 \cdot (1 - p)}{p}

Equivalently, the tumor's deviation from 0.5 is the observed deviation
amplified by :math:`1/p`: when purity is low, small mis-modeled deviations
in the observed BAF are amplified into large deviations in the rescaled
output.

When the input matches the model assumptions, the rescaled tumor BAF lies in
[0, 1]. If the rescaled value falls outside that interval, the values are
clamped to [0, 1] and a warning is logged. See :ref:`baf-troubleshooting`.


Allele-specific copy number and LOH
-----------------------------------

The :ref:`call` command pairs the rescaled BAF with the integer total copy
number to split each segment into its two allelic copy numbers, output as
columns ``cn1`` (major) and ``cn2`` (minor) in the ``.cns`` file. The
calculation follows `PSCBS <http://bioinformatics.oxfordjournals.org/content/27/15/2038.short>`_:
the total copy number is multiplied by the upper-half BAF and rounded to the
nearest integer, with the constraint :math:`\text{cn1} \geq \text{cn2}` and
:math:`\text{cn1} + \text{cn2} = \text{cn}`.

Allelic imbalance, including copy-number-neutral loss of heterozygosity
(LOH), is apparent when ``cn1`` and ``cn2`` differ. Specifically:

- ``cn1 == cn2``: balanced segment (e.g. ``2/2`` for diploid neutral, ``3/3``
  for a balanced 6-copy gain).
- ``cn1 != cn2``: allelic imbalance (e.g. ``2/1`` for a single-copy loss).
- ``cn2 == 0`` with ``cn > 0``: complete loss of heterozygosity for that
  segment (e.g. ``2/0`` for copy-number-neutral LOH, ``1/0`` for hemizygous
  loss).

If the segment had no overlapping heterozygous SNPs, the ``baf``, ``cn1``,
and ``cn2`` columns are written as missing values (NaN).


.. _baf-troubleshooting:

Troubleshooting
---------------

The two most common BAF problems are summarized below.

Negative or out-of-range BAF in ``.cns`` output
```````````````````````````````````````````````

Symptom: ``cnvkit call --purity`` produces a ``.cns`` file with ``baf``
values outside [0, 1] (older CNVkit versions), or you see a log message
like::

    WARNING: 17 segment(s) had tumor BAF outside [0, 1] after purity
    rescaling (purity=0.37); values clamped. The purity estimate may be
    too low, or the input VCF may contain somatic variants.

This means observed BAFs were too far from 0.5 for the rescaling model
(see :ref:`baf-rescaling`). Likely causes:

- **The VCF contains somatic-only variants** (most common). Re-run with a
  VCF that includes germline heterozygous SNPs. See :ref:`baf-vcf-prep`.
- **Sample IDs not specified.** With ``-i``/``-n`` missing in a paired
  tumor/normal VCF, CNVkit cannot apply genotype-based somatic filtering.
  Set ``--sample-id`` and ``--normal-id`` (or add a PEDIGREE header).
- **Purity estimate is too low.** If the purity passed to ``--purity`` is
  below the actual tumor cell fraction, observed deviations from 0.5 get
  amplified beyond the [0, 1] range during rescaling.

Out-of-range values are clamped to [0, 1] on output, but the underlying
input issue should be investigated.

"No heterozygous variants" or "Median allele frequency far from 0.5"
````````````````````````````````````````````````````````````````````

Symptom: a log message at load time, e.g.::

    WARNING: No heterozygous variants remain after filtering.

or::

    WARNING: Median allele frequency 0.93 is far from the 0.5 expected for
    heterozygous germline SNPs.

Both indicate the loaded variants don't look like germline heterozygous
SNPs. Re-check the VCF preparation (see :ref:`baf-vcf-prep`) and sample
identification.
