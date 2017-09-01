Plots and graphics
==================

.. _scatter:

scatter
-------

Plot bin-level log2 coverages and segmentation calls together.  Without any
further arguments, this plots the genome-wide copy number in a form familiar to
those who have used array CGH.

::

    cnvkit.py scatter Sample.cnr -s Sample.cns
    # Shell shorthand
    cnvkit.py scatter -s TR_95_T.cn{s,r}

.. image:: TR_95_T-scatter.png


The options ``--chromosome`` and ``--gene`` (or their single-letter equivalents)
focus the plot on the specified region::

    cnvkit.py scatter -s Sample.cn{s,r} -c chr7
    cnvkit.py scatter -s Sample.cn{s,r} -c chr7:140434347-140624540
    cnvkit.py scatter -s Sample.cn{s,r} -g BRAF

In the latter two cases, the genes in the specified region or with the specified
names will be highlighted and labeled in the plot. The ``--width`` (``-w``)
argument determines the size of the chromosomal regions to show flanking the
selected region. Note that only targeted genes can be highlighted and labeled;
genes that are not included in the list of targets are not labeled in the .cnn
or .cnr files and are therefore invisible to CNVkit.

The arguments ``-c`` and ``-g`` can be combined to e.g. highlight
specific genes in a larger context::

    # Show a chromosome arm, highlight one gene
    cnvkit.py scatter -s Sample.cn{s,r} -c chr5:100-50000000 -g TERT
    # Show the whole chromosome, highlight two genes
    cnvkit.py scatter -s Sample.cn{s,r} -c chr7 -g BRAF,MET
    # Highlight two genes in a specified range
    cnvkit.py scatter -s TR_95_T.cn{s,r} -c chr12:50000000-80000000 -g CDK4,MDM2

.. image:: TR_95_T-CDK4-MDM2-scatter.png

When a chromosomal region is plotted with CNVkit's "scatter" command , the size
of the plotted datapoints is proportional to the weight of each point used in
segmentation -- a relatively small point indicates a less reliable bin.
Therefore, if you see a cluster of smaller points in a short segment (or where
you think there ought to be a segment, but there isn't one), then you can cast
some doubt on the copy number call in that region. The dispersion of points
around the segmentation line also visually indicates the level of noise or
uncertainty.

To create multiple region-specific plots at once, the regions of interest can be
listed in a separate file and passed to the ``scatter`` command with the
``-l``/``--range-list`` option. This is equivalent to creating the plots
separately with the ``-c`` option and then combining the plots into a single
multi-page PDF.

The bin-level log2 ratios or coverages can also be plotted without segmentation
calls::

    cnvkit.py scatter Sample.cnr

This can be useful for viewing the raw, un-corrected coverage depths when
deciding which samples to use to build a profile, or simply to see the coverages
without being helped/biased by the called segments.

The ``--trend`` option (``-t``) adds a smoothed trendline to the plot. This is
fairly superfluous if a valid segment file is given, but could be helpful if the
CBS dependency is not available, or if you're skeptical of the segmentation in a
region.

SNV b-allele frequencies
````````````````````````

Loss of heterozygosity (LOH) can be viewed alongside copy number by passing
variants as a VCF file with the ``-v`` option. Heterozygous SNP allelic
frequencies are shown in a subplot below the CNV scatter plot.

::

    cnvkit.py scatter Sample.cnr -s Sample.cns -v Sample.vcf

If only the VCF file is given by itself, just plot the allelic frequencies::

    cnvkit.py scatter -v Sample.vcf

Given segments, show the mean b-allele frequency values above and below 0.5 of
SNVs falling within each segment. Divergence from 0.5 indicates LOH in the tumor
sample.

::

    cnvkit.py scatter -s Sample.cns -v Sample.vcf -i TumorID -n NormalID

Regions with LOH are reflected in heterozygous germline SNPs in the tumor sample
with allele frequencies shifted away from the expected 0.5 value.
Given a VCF with only the tumor sample called, it is difficult to focus on just
the informative SNPs because it's not known which SNVs are present and
heterozygous in normal, germline cells.
Better results can be had by giving CNVkit more information:

- Call somatic mutations using paired tumor and normal samples.
  In the VCF, the somatic variants should be flagged in the INFO column with the
  string "SOMATIC". (MuTect does this automatically.) Then CNVkit will skip
  these for plotting.
- Add a "PEDIGREE" tag to the VCF header, listing the tumor sample as "Derived"
  and the normal as "Original". (MuTect doesn't do this, but it does add a
  nonstandard GATK header that CNVkit can extract the same information from.)
- In lieu of a PEDIGREE tag, tell CNVkit which sample IDs are the tumor and normal using the
  ``-i`` and ``-n`` options, respectively.
- If no paired normal sample is available, you can still filter for likely
  informative SNPs by intersecting your tumor VCF with a set of known SNPs such
  as 1000 Genomes, ESP6500, or ExAC.
  Drop the private SNVs that don't appear in these databases to create a VCF
  more amenable to LOH detection.


.. _diagram:

diagram
-------

Draw copy number (either individual bins (.cnn, .cnr) or segments (.cns)) on
chromosomes as an ideogram. If both the bin-level log2 ratios and segmentation
calls are given, show them side-by-side on each chromosome (segments on the left
side, bins on the right side).

::

    cnvkit.py diagram Sample.cnr
    cnvkit.py diagram -s Sample.cns
    cnvkit.py diagram -s Sample.cns Sample.cnr

If bin-level log2 ratios are provided (.cnr), genes with log2 ratio values
beyond a fixed threshold will be labeled on the plot.
This plot style works best with target panels of a few hundred genes at most;
with whole-exome sequencing there are often so many genes affected by CNAs that
the individual gene labels become difficult to read.

.. image:: TR_95_T-diagram.png

If only segments are provided (``-s``), gene labels are not shown.  This plot is
then equivalent to the ``heatmap`` command, which effectively summarizes the
segmented values from many samples.

By default, the sex chromosomes X and Y are colorized relative to the expected
ploidy, i.e. for female samples analyzed with a male reference, while the X
chromosome has a copy ratio near +1.0 in the input .cnr and .cns files, in the
output diagram it will be shown as neutral copy number (white or faint colors)
rather than a gain (red), because the diploid X is expected. The sample sex can
be specified with the ``-x``/``--sample-sex`` option, or will otherwise be
guessed automatically (see :doc:`sex`). This correction is done by default, but
can be disabled with the option ``--no-shift-xy``.


.. _heatmap:

heatmap
-------

Draw copy number (either bins (.cnn, .cnr) or segments (.cns)) for multiple
samples as a heatmap.

To get an overview of the larger-scale CNVs in a cohort, use the
"heatmap" command on all .cns files::

    cnvkit.py heatmap *.cns

.. image:: heatmap-tr-nod.png

The color range can be subtly rescaled with the ``-d`` option to de-emphasize
low-amplitude segments, which are likely spurious CNAs::

    cnvkit.py heatmap *.cns -d

.. image:: heatmap-tr.png

A heatmap can also be drawn from bin-level log2 coverages or copy ratios (.cnn,
.cnr), but this will be extremely slow at the genome-wide level.
Consider doing this with a smaller number of samples and only for one chromosome
or chromosomal region at a time, using the ``-c`` option::

    cnvkit.py heatmap TR_9*T.cnr -c chr12  # Slow!
    cnvkit.py heatmap TR_9*T.cnr -c chr7:125000000-145000000

.. image:: heatmap-tr-chr12.png

If an output file name is not specified with the ``-o`` option, an interactive
matplotlib window will open, allowing you to select smaller regions, zoom in,
and save the image as a PDF or PNG file.

The samples are shown in the order there's given on the command line.
If you use "\*.cns" then the filenames might always be fetched
alphabetically (depending on your operating system), but if you type
them out in the order you like, it should keep that order. You can use
the Unix shell to pull the names out of a file on the fly, e.g.::

    cnvkit.py heatmap `cat filenames.txt`


As with :ref:`diagram`, the sex chromosomes X and Y are colorized relative to
the expected ploidy, based on the sample and reference sex (see :doc:`sex`).
This correction can be disabled with the option ``--no-shift-xy``.


.. _plotcustom:

Customizing plots
-----------------

The plots generated with the :ref:`scatter` and :ref:`heatmap` commands use the
Python plotting library matplotlib.

To quickly adjust the displayed area of the genome in a plot, run either
plotting command without the `-o` option to generate an interactive plot in a
new window. You can then resize that plot up to the full size of your screen,
use the plot window's selection mode to select a smaller area of the genome, and
use the plot window's save button to save the plot in your preferred format.

You can customize font sizes and other aspects of the plots by `configuring
matplotlib <http://matplotlib.org/users/customizing.html>`_.
If you're running CNVkit on the command line and not using it as a Python
library, then you can just create a file in your home directory (or the same
directory as ``cnvkit.py``) called ``.matplotlibrc``.  For example, to shrink
the font size of the x- and y-axis labels, put this line in the configuration
file::

    axes.labelsize      : small

For more control, in the Python intepreter (or a script, or a Jupyter notebook),
import the :doc:`cnvlib` module and call the `do_scatter` or `do_heatmap`
function to create a plot. Then you can use matplotlib.pyplot to get the current
axis and modify the plot elements, change font sizes, or anything else you
like::

    from glob import glob
    from matplotlib import pyplot as plt
    import cnvlib

    segments = map(cnvlib.read, glob("*.cns"))
    ax = cnvlib.do_heatmap(segments)
    ax.set_title("All my samples")
    plt.rcParams["font.size"] = 9.0
    plt.show()
