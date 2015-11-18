Compatibility and other I/O
===========================

.. _import-picard:

import-picard
-------------

Convert Picard CalculateHsMetrics per-target coverage files (.csv) to the
CNVkit .cnn format::

    cnvkit.py import-picard *.hsmetrics.targetcoverages.csv *.hsmetrics.antitargetcoverages.csv
    cnvkit.py import-picard picard-hsmetrics/ -d cnvkit-from-picard/

You can use `Picard tools <http://broadinstitute.github.io/picard/>`_ to perform
the bin read depth and GC calculations that CNVkit normally performs with the
:ref:`coverage` and :ref:`reference` commands, if need be.

Procedure:

1. Use the :ref:`target` and :ref:`antitarget` commands to generate the
   "targets.bed" and "antitargets.bed" files.
2. Convert those BED files to Picard's "interval list" format by adding the BAM
   header to the top of the BED file and rearranging the columns -- see the
   Picard command `BedToIntervalList
   <http://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList>`_.
3. Run Picard `CalculateHsMetrics
   <http://broadinstitute.github.io/picard/command-line-overview.html#CalculateHsMetrics>`_
   on each of your normal/control BAM files with the "targets" and "antitargets"
   interval lists (separately), your reference genome, and the
   "PER_TARGET_COVERAGE" option.
4. Use :ref:`import-picard` to convert all of the PER_TARGET_COVERAGE files to
   CNVkit's .cnn format.
5. Use :ref:`reference` to build a CNVkit reference from those .cnn files. It
   will retain the GC values Picard calculated; you don't need to provide the
   reference genome sequence again to get GC (but you if you do, it will also
   calculate the RepeatMaster fraction values)
6. Use :ref:`batch` with the ``-r``/``--reference`` option to process the rest
   of your test samples.


.. _import-seg:

import-seg
----------

Convert a file in the `SEG format <https://www.broadinstitute.org/igv/SEG>`_
(e.g. the output of standard CBS or the GenePattern server) into one or more
CNVkit .cns files.

The chromosomes in a SEG file may have been converted from chromosome names to
integer IDs. Options in ``import-seg`` can help recover the original names.

* To add a "chr" prefix, use "-p chr".
* To convert chromosome indices 23, 24 and 25 to the names "X", "Y" and "M" (a
  common convention), use "-c human".
* To use an arbitrary mapping of indices to chromosome names, use a
  comma-separated "key:value" string. For example, the human convention would
  be: "-c 23:X,24:Y,25:M".


.. _import-theta:

import-theta
------------

Convert the ".results" output of `THetA2
<http://compbio.cs.brown.edu/projects/theta/>`_ to one or more CNVkit .cns files
representing subclones with integer absolute copy number in each segment.


.. _export:

export
------

Convert copy number ratio tables (.cnr files) or segments (.cns) to
another format.

bed
```

Segments can be exported to BED format to support a variety of other uses, such
as viewing in a genome browser.  The log2 ratio value of each segment is
converted and rounded to an integer value, as required by the BED format. To get
accurate copy number values, see the :ref:`call` command.

::

    # Estimate integer copy number of each segment
    cnvkit.py call Sample.cns -y -o Sample.call.cns
    # Show estimated integer copy number of all regions
    cnvkit.py export bed Sample.call.cns --show all -y -o Sample.bed

The same format can also specify CNV regions to the FreeBayes variant caller
with FreeBayes's ``--cnv-map`` option::

    # Show only CNV regions
    cnvkit.py export bed Sample.call.cns -o all-samples.cnv-map.bed

By default only regions with copy number different from the given ploidy
(default 2) are output. (Notice what this means for allosomes.)
To output all segments, use the ``--show all`` option.

vcf
```

Convert segments, ideally already adjusted by the :ref:`call` command, to
a :ref:`vcfformat` file. Copy ratios are converted to absolute integers, as with
BED export, and VCF records are created for the segments where the copy number
is different from the expected ploidy (e.g. 2 on autosomes, 1 on haploid sex
chromosomes, depending on sample gender).

Gender can be specified with the ``-g``/``--gender`` option, or will be guessed
automatically. If a male reference is used, use ``-y``/``--male-reference`` to
say so. Note that these are different: If a female sample is run with a male
reference, segments on chromosome X with log2-ratio +1 will be skipped, because
that's the expected copy number, while an X-chromosome segment with log2-ratio 0
will be printed as a hemizygous loss.

::

    cnvkit.py export vcf Sample.cns -y -g female -i "SampleID" -o Sample.cnv.vcf

cdt, jtv
````````

A collection of probe-level copy ratio files (``*.cnr``) can be exported to Java
TreeView via the standard CDT format or a plain text table::

    cnvkit.py export jtv *.cnr -o Samples-JTV.txt
    cnvkit.py export cdt *.cnr -o Samples.cdt

seg
```

Similarly, the segmentation files for multiple samples (``*.cns``) can be
exported to the standard SEG format to be loaded in the Integrative Genomic
Viewer (IGV)::

    cnvkit.py export seg *.cns -o Samples.seg

nexus-basic
```````````

The format ``nexus-basic`` can be loaded directly by the commercial program
Biodiscovery Nexus Copy Number, specifying the "basic" input format in that
program. This allows viewing CNVkit data as if it were from array CGH.

This is a tabular format very similar to .cnr files, with the columns:

#. chromosome
#. start
#. end
#. log2


nexus-ogt
`````````

The format ``nexus-ogt`` can be loaded directly by the commercial program
Biodiscovery Nexus Copy Number, specifying the "Custom-OGT" input format in that
program. This allows viewing CNVkit data as if it were from a SNP array.

This is a tabular format similar to .cnr files, but with B-allele frequencies
(BAFs) extracted from a corresponding VCF file. The format's columns are (with
.cnr equivalents):

#. "Chromosome" (chromosome)
#. "Position" (start)
#. "Position" (end)
#. "Log R Ratio" (log2)
#. "B-Allele Frequency" (from VCF)

The positions of each heterozygous variant record in the given VCF are matched
to bins in the given .cnr file, and the variant allele frequencies are extracted
and assigned to the matching bins.

- If a bin contains no variants, the BAF field is left blank
- If a bin contains multiple variants, the BAFs of those variants are "mirrored"
  to be all above .5 (e.g. BAF of .3 becomes .7), then the median is taken as
  the bin-wide BAF.

