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

Convert a file in the SEG format (e.g. the output of standard CBS or the
GenePattern server) into one or more CNVkit .cns files.

The chromosomes in a SEG file may have been converted from chromosome names to
integer IDs. Options in ``import-seg`` can help recover the original names.

* To add a "chr" prefix, use "-p chr".
* To convert chromosome indices 23, 24 and 25 to the names "X", "Y" and "M" (a
  common convention), use "-c human".
* To use an arbitrary mapping of indices to chromosome names, use a
  comma-separated "key:value" string. For example, the human convention would
  be: "-c 23:X,24:Y,25:M".


.. _export:

export
------

Convert copy number ratio tables (.cnr files) to another format.

bed
```

The segmented output from multiple samples (``*.cns``) can be exported to BED
format to support a variety of other uses, such as viewing in a genome browser.
The log2 ratio value of each segment is converted and rounded to an integer
value, as required by the BED format. To get accurate copy number values, see
the :ref:`call` command.

::

    # Estimate integer copy number of each segment
    cnvkit.py call Sample.cns -y -o Sample.call.cns
    # Show estimated integer copy number of all regions
    cnvkit.py export bed Sample.call.cns --show-neutral -y -o Sample.bed

The same format can also specify CNV regions to the FreeBayes variant caller
with FreeBayes's ``--cnv-map`` option::

    # Show only CNV regions
    cnvkit.py export bed *.call.cns -o all-samples.cnv-map.bed

Copy-number-neutral regions are not shown in the output by default, but can be
included with the ``--show-neutral`` option.

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
program.

