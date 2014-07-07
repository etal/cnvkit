Compatibility and other I/O
===========================

import-picard
-------------

Convert Picard CalculateHsMetrics coverage files (.csv) to the CNVkit .cnn
format.

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

export
------

Convert copy number ratio tables (.cnr files) to another format.

A collection of probe-level copy ratio files (``*.cnr``) can be exported to Java
TreeView via the standard CDT format or a plain text table::

    cnvkit.py export jtv *.cnr -o Samples-JTV.txt
    cnvkit.py export cdt *.cnr -o Samples.cdt

Similarly, the segmentation files for multiple samples (``*.cns``) can be
exported to the standard SEG format to be loaded in the Integrative Genomic
Viewer (IGV)::

    cnvkit.py export seg *.cns -o Samples.seg

Also note that the individual ``.cnr`` and ``.cnn`` files can be loaded directly
by the commercial program Biodiscovery Nexus Copy Number, specifying the "basic"
input format.

