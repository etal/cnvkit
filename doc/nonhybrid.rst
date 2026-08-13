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
Since the input does not contain useful per-target gene labels, a gene
annotation database is recommended to label genes in the outputs::

    cnvkit.py batch Sample1.bam Sample2.bam -n Control1.bam Control2.bam \
            -m wgs -f hg38.fasta --annotate data/refFlat_hg38.txt

To speed up and/or improve the accuracy of WGS analyses, try any or all of the
following:

- Instead of analyzing the whole genome, use the "target" BED file
  to limit the analysis to just the genic regions. You can get such a BED file
  from the `UCSC Genome Browser <https://genome.ucsc.edu/cgi-bin/hgTables>`_, for
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


Distributing a cohort across a cluster
``````````````````````````````````````

On a single machine, run the whole cohort in one ``batch`` command: pass every
tumor BAM at once and let ``-p``/``--processes`` divide the available CPUs
among them. Splitting the cohort into separate commands is worthwhile only to
reach more than one machine, since no single invocation can span nodes. To do
that, build the pooled reference once, then analyze each tumor sample as an
independent job against that reference.

The first job builds the reference from the normal samples alone. With no tumor
BAMs on the command line, ``batch`` stops after writing the reference::

    cnvkit.py batch -n Control*.bam -m wgs -f hg38.fasta \
            --annotate data/refFlat_hg38.txt \
            --output-reference cohort/ref-wgs.cnn -d cohort/ -p 8

Here ``-f``/``--fasta`` supplies both the sequencing-accessible regions, which
are computed on the fly, and the GC and repeat-masking corrections; give
``-g``/``--access`` as well if you have already run the :ref:`access` command.

Each subsequent job analyzes one tumor sample against the finished reference::

    cnvkit.py batch Sample1.bam -m wgs -r cohort/ref-wgs.cnn \
            -d cohort/Sample1 -p 8

Submit one such job per tumor sample. The jobs are independent -- each reads
the reference and writes only outputs named after its own sample, the full set
being listed under :ref:`batch` -- so they can be run in any order and on
separate nodes. Pass ``-m wgs`` to these jobs as well: the sequencing method
selects the WGS bias corrections and segmentation threshold for the sample, and
is not recorded in the reference file. The options that construct a reference,
including ``-n``, ``-t``, ``--annotate`` and ``--target-avg-size``, are rejected
when ``-r`` is given.

Two practical points:

- Give each sample its own output directory, which ``batch`` creates if it does
  not exist. Besides keeping the per-sample outputs apart, this avoids
  concurrent jobs rewriting one shared file: ``batch`` extracts the reference's
  bin coordinates into temporary target and antitarget BED files in the output
  directory, named after the reference because their contents depend only on
  it, and a job re-writing them while another job is reading them can be seen
  mid-write.
- ``-p``/``--processes`` parallelizes work within a single job, across the CPUs
  of the node it lands on; it does not distribute work between nodes. ``-p 0``
  uses every CPU the job is permitted to use, honoring the CPU affinity or
  cgroup limits imposed by the scheduler.

``batch`` does not checkpoint. It keeps no run state between invocations and
does not skip a step whose output file already exists, so re-running a sample
recomputes that sample's coverage, bin-level copy ratios and segments from the
BAM. The unit of loss is therefore a single sample: when one job fails, re-run
that job alone. Samples that already finished are unaffected, since nothing
outside their own job reads or rewrites their outputs, and the reference does
not need to be rebuilt, since it is an input to the per-sample jobs rather than
an output of them. The restarted sample begins again at coverage calculation,
which for WGS is the longest step.

.. note::
    **Index every BAM before submitting the jobs.** CNVkit indexes a BAM
    automatically when the index is missing or older than the BAM, but it
    writes the index next to the BAM file: if that directory is not writable,
    the run stops with a samtools error reporting that it failed to create or
    write the index. Indexing the cohort in advance with ``samtools index``
    avoids this, and keeps the indexing time out of the cluster jobs. See
    :ref:`coverage` for the BAM preparation requirements in full.


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
    # For each sample
    cnvkit.py coverage Sample.bam targets.split.bed -p 0 -o Sample.targetcoverage.cnn
    cnvkit.py reference *.targetcoverage.cnn --no-edge -o ref-tas.cnn
    cnvkit.py fix Sample.targetcoverage.cnn -r ref-tas.cnn --no-edge

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
    <https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/5774-removing-duplicates-is-it-really-necessary>`_
    has details and background).

    In targeted amplicon sequencing, all of the amplified reads are in fact PCR
    duplicates by design. By marking and thus omitting these reads, the
    remaining coverage will be low, as if no amplification were performed.
