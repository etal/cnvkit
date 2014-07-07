Plots and graphics
==================

scatter
-------

Plot probe log2 coverages and segmentation calls together.

::

    cnvkit.py scatter Sample.cnr -s Sample.cns

The options ``--gene``, ``--chromosome`` or ``--range`` (or their single-letter
equivalents) focus the plot on the specified region::

    cnvkit.py scatter Sample.cnr -s Sample.cns -r chr7
    cnvkit.py scatter Sample.cnr -s Sample.cns -r BRAF
    cnvkit.py scatter Sample.cnr -s Sample.cns -r chr7:140434347-140624540

In the latter two cases, the ``--width`` (``-w``) argument determines the size
of the chromosomal regions to show flanking the selected region.

Loss of heterozygosity (LOH) can be viewed alongside copy number by passing
variants as a VCF file with the ``-v`` option. Heterozygous SNP allelic
frequencies are shown in a subplot below the CNV scatter plot. (Also see the
``loh`` command, below.)

::

    cnvkit.py scatter Sample.cnr -s Sample.cns -v Sample.vcf

The probe copy number values can also be plotted without segmentation calls::

    cnvkit.py scatter Sample.cnr

This can be useful if CBS is unavailable, or for viewing the raw, un-corrected
coverages when deciding which samples to use to build a profile, or simply to
see the coverages without being helped/biased by the called segments.

The ``--trend`` option (``-t``) adds a smoothed trendline to the plot. This is
fairly superfluous if a valid segment file is given, but could be helpful if CBS
is not available, or if you're skeptical of the segmentation in a region.

loh
---

Plot allelic frequencies at each variant position in a VCF file. Divergence from
0.5 indicates loss of heterozygosity (LOH) in a tumor sample.

::

    cnvkit.py loh Sample.vcf

diagram
-------

Draw copy number (either raw probes (.cnn, .cnr) or segments (.cns)) on
chromosomes as a diagram. If both the raw probes and segmentation calls are
given, show them side-by-side on each chromosome (segments on the left side,
probes on the right side).

::

    cnvkit.py diagram Sample.cnr
    cnvkit.py diagram -s Sample.cns
    cnvkit.py diagram -s Sample.cns Sample.cnr

heatmap
-------

Draw copy number (either raw probes (.cnn, .cnr) or segments (.cns)) for
multiple samples as a heatmap.

The segmentation calls alone will render much faster, and will probably be more
useful to look at.

::

    cnvkit.py heatmap *.cns
    cnvkit.py heatmap *.cnr  # Slow!

