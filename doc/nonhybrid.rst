Whole-genome sequencing and targeted amplicon capture
=====================================================

CNVkit is primarily designed for use on **hybrid capture** sequencing data,
where off-target reads are present and can be used improve copy number
estimates. However, CNVkit can also be used on **whole-genome sequencing** (WGS)
and **targeted amplicon sequencing** (TAS) datasets by using alternative
command-line options.

The :ref:`batch` command supports these workflows through the
``-m``/``--method`` option.


.. _wgs:

Whole-Genome Sequencing (WGS)
-----------------------------

CNVkit treats WGS data as a capture of all of the genome's sequencing-accessible
regions, with no off-target regions.

The ``batch --method wgs`` option uses the given reference genome's
sequencing-accessible regions ("access" BED) as the "targets" -- these will be
calculated on the fly if not provided. No "antitarget" regions are used.
Since the input does not contain useful per-target gene labels, a  gene
annotation database is required and used to label genes in the outputs::

    cnvkit.py batch Sample1.bam Sample2.bam -n Control1.bam Control2.bam \
            -m wgs -f hg19.fasta --annotate refFlat.txt

To speed up and/or improve the accuracy of WGS analyses, try any or all of the
following:

- Instead of analyzing the whole genome, use the "target" BED file
  to limit the analysis to just the genic regions. You can get such a BED file
  from the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables), for
  example.
- Increase the "target" average bin size (``--target-avg-size``), e.g. to at
  least 1000 bases for 30x coverage, or proportionally more for lower-coverage
  sequencing.
- Specify a smaller p-value threshold (``segment -t``). For the CBS method,
  ``1e-6`` may work well. Or, try the ``hmm`` segmentation method.
- Use the ``-p/--processes`` option in the :ref:`batch`, :ref:`coverage` and
  :ref:`segment` commands to ensure all available CPUs are used.
- Ensure you are using the most recent version of CNVkit. Each release includes
  some performance improvements.
- Turn off the "edge" bias correction in the :ref:`reference` and :ref:`fix`
  commands (`--no-edge`).

The ``batch -m wgs`` option does all of these except the first automatically.


.. _tas:

Targeted Amplicon Sequencing (TAS)
----------------------------------

When amplicon sequencing is used as a targeted capture method, no off-target
reads are sequenced. While this limits the copy number information available in
the sequencing data versus hybrid capture, CNVkit can analyze TAS data using
only on-target coverages and excluding all off-target regions from the analysis.

The ``batch -m amplicon`` option uses the given targets to infer coverage,
ignoring off-target regions::

    cnvkit.py batch -m amplicon -t targets.bed *.bam

Equivalently::

    cnvkit.py target targets.bed --split -o targets.split.bed
    # Create a blank file to substitute for antitargets
    touch MT
    # For each sample
    cnvkit.py coverage Sample.bam targets.split.bed -p 0 -o Sample.targetcoverage.cnn
    cnvkit.py reference *.targetcoverage.cnn --no-edge -o ref-tas.cnn
    cnvkit.py fix Sample.targetcoverage.cnn MT ref-tas.cnn --no-edge

This approach does not collect any copy number information between targeted
regions, so it should only be used if you have in fact prepared your samples
with a targeted amplicon sequencing protocol. It also does not attempt to
further normalize each amplicon at the gene level, though this may be addressed
in a future version of CNVkit.

.. note::
    **Do not mark duplicates** in the BAM files for samples sequenced by this
    method.

    Picard MarkDuplicates, samtools rmdup, *et al.* are designed to flag
    possible PCR duplicates (originally for WGS datasets, but also useful for
    hybrid capture). Variant callers like GATK and CNVkit will ignore those
    reads in their internal calculations, considering these reads to be
    non-independent measurements. (`This SeqAnswers thread
    <http://seqanswers.com/forums/showthread.php?t=6854>`_ has details and
    background).

    In targeted amplicon sequencing, all of the amplified reads are in fact PCR
    duplicates by design. By marking and thus omitting these reads, the
    remaining coverage will be low, as if no amplification were performed.
