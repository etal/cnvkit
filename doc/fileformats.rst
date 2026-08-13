File formats
============

We've tried to use standard file formats where possible in CNVkit. However, in a
few cases we have needed to extend the standard BED format to accommodate
additional information.

All of the non-standard file formats used by CNVkit are tab-separated plain text
and can be loaded in a spreadsheet program, R or other statistical analysis
software for manual analysis, if desired.

Every interval CNVkit reads must end at or after it starts. In a region file
written by another tool, a reversed pair is taken to span the wider range and
the repair is logged. For BED that is what the writer meant, since some
primer-design tools record reverse-direction amplicons that way; for an
annotation file carrying its own strand column, such as GFF or refFlat, it is
instead the only usable interval left in a damaged row. In
CNVkit's own tables below, the ``.cnn``, ``.cnr``, ``.cns`` and SEG formats, a
reversed row is refused: those coordinates came from a binner or a
segmenter, so a reversed one means the file was damaged after it was written,
and repairing it would invent a span no segmenter produced. Versions before
0.9.14 could write such a row themselves when squashing segments whose bins
were out of order (#677), so regenerating an old file with a current version
is usually enough.

.. _bedformat:

BED and GATK/Picard Interval List
---------------------------------

- UCSC Genome Browser's `BED definition and FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_
- GATK's `Interval List description
  <https://www.broadinstitute.org/gatk/guide/article?id=1204>`_ and `FAQ
  <https://www.broadinstitute.org/gatk/guide/article?id=1319>`_

Note that BED genomic coordinates are 0-indexed, like C or Python code -- for
example, the first nucleotide of a 1000-basepair sequence has position 0, the
last nucleotide has position 999, and the entire region is indicated by the
range 0-1000.

GATK and Picard interval list coordinates are 1-indexed, like R or Matlab code.
In the same example, the first nucleotide of a 1000-basepair sequence has
position 1, the last nucleotide has position 1000, and the entire region is
indicated by the range 1-1000. These files usually have the extension
``.interval_list``.

In GATK4, the term "interval list" also refers to samtools-style genomic
coordinate specifications of the form *chromosome:start-end*, e.g.
``chr1:1-1000``. As with Picard and older GATK style interval lists, the
coordinates are 1-indexed. When used with GATK4, these files usually have the
extension ``.list`` or ``.interval``.

CNVkit will load these files by automatically determining the specific format
based on the file contents, not the filename extension.


.. _gffformat:

GFF
---

CNVkit can read `GFF3 <http://gmod.org/wiki/GFF3>`_, `GFF2
<http://gmod.org/wiki/GFF2>`_, and `GTF <http://mblab.wustl.edu/GTF2.html>`_
files as input in most commands where UCSC BED works. These formats all have these 9 columns:

    - seqid/reference/seqname/chromosome
    - source
    - type/method/feature
    - start: in 1-based integer coordinates
    - end: in 1-based integer coordinates
    - score: float or '.' (for NA)
    - strand: [+-.?]
    - phase/frame: [012.]
    - attribute/group: string

The difference between the formats is in column 9; since CNVkit only attempts
to extract the gene name (at most) from this column, the other features are
ignored and the formats are effectively the same for the purpose of labeling
genomic regions. These docs therefore refer to these formats collectively as GFF.

To extract gene names, CNVkit's GFF reader checks for these tags in column 9 in order:

    - ``Name``
    - ``gene_id``
    - ``gene_name``
    - ``gene``

It will take the value of the first match and use it as the "gene" in the
internal data structures and in other output file formats.


 .. _segformat:

SEG
---

The SEG format is the `tabular output
<https://software.broadinstitute.org/software/igv/SEG>`_ of DNAcopy, the
reference implementation of Circular Binary Segmentation (CBS). It is a
tab-separated table with the following 5 or 6 columns:

    - ``ID`` -- sample name
    - ``chrom`` -- chromosome name or ID
    - ``loc.start`` -- segment's genomic start position, 1-indexed
    - ``loc.end`` -- segment end position
    - ``num.mark`` -- (optional) number of probes or bins covered by the segment
    - ``seg.mean`` -- segment mean value, usually in log2 scale

The column names in the first line are not enforced, and can vary across
implementations.

SEG files can be used with a number of other programs that operate on segmented
log2 copy ratios -- including GISTIC 2.0, IGV, the GenePattern server, and many
R packages.

To convert CNVkit's .cns files to SEG, use the command :ref:`export` ``seg``,
and to convert SEG files produced outside of CNVkit into CNVkit's own segmented
format (.cns), use :ref:`import-seg`.


.. _vcfformat:

VCF
---

See the `VCF specifications <https://github.com/samtools/hts-specs>`_.

CNVkit currently uses VCF files in two ways:

- To extract single-nucleotide variant (SNV) allele frequencies, which can be
  plotted in the :ref:`scatter` command, used to assign allele-specific copy
  number in the :ref:`call` command, or exported along with bin-level copy
  ratios to the "nexus-ogt" format.
- To :ref:`export` CNVs, describing/encoding each CNV segment as a structural
  variant (SV).

For the former -- investigating allelic imbalance and loss of heterozygosity
(LOH) -- the VCF must contain heterozygous **germline** SNVs and, ideally,
matched tumor/normal samples with the tumor identified via the PEDIGREE
header tag or ``-i``/``-n`` command-line options. See :doc:`baf` for full
guidance on VCF preparation, sample identification, and troubleshooting.

CNVkit reads each sample's zygosity from the ``GT`` genotype field when it is
present: a call with both a reference and an alternate allele (e.g. ``0/1``) is
heterozygous, ``0/0`` is homozygous reference, and ``1/1`` is homozygous
alternate. A no-call allele (``.``) is treated as missing rather than as a
distinct allele, so a complete no-call (``./.``) and a genotype with no ``GT``
field at all (e.g. Strelka somatic output) fall back to inferring zygosity from
the alternate-allele frequency. A partial call keeps its called allele: ``0/.``
(no alternate observed) is homozygous reference, while ``1/.`` retains the
alternate allele and is treated as heterozygous.


.. _cnxformat:

Target and antitarget bin-level coverages (.cnn)
------------------------------------------------

CNVkit saves its information in a tabular format similar to BED, but with
additional columns.  Each row in the file indicates an on-target or off-target
(a.k.a. "antitarget") bin. Genomic coordinates are 0-indexed, like BED.
Column names are shown as the first line of the file.

In the output of the :ref:`coverage` command, the columns are:

* Chromosome or reference sequence name (``chromosome``)
* Start position (``start``)
* End position (``end``)
* Gene name (``gene``)
* Log2 mean coverage depth (``log2``)
* Absolute-scale mean coverage depth (``depth``)

Essentially the same tabular file format is used for coverages (.cnn), ratios
(.cnr) and segments (.cns) emitted by CNVkit.


Copy number reference profile (.cnn)
------------------------------------

In addition to the columns present in the "target" and "antitarget" .cnn files,
the reference .cnn file has the columns:

* GC content of the sequence region (``gc``)
* RepeatMasker-masked proportion of the sequence region (``rmask``)
* Statistical spread or dispersion (``spread``)

The **log2** coverage depth is the robust average of coverage depths,
excluding extreme outliers, observed at the corresponding bin in each the sample
.cnn files used to construct the :ref:`reference`. The **spread** is a similarly
robust estimate of the standard deviation of normalized log2 coverages in the
bin. The **depth** column is the robust average of absolute-scale coverage
depths from the input .cnn files, but without any bias corrections.

To manually review potentially problematic targets in the built reference, you
can sort the file by the **spread** column; bins with higher values are the
noisy ones.

It is important to keep the copy number reference file consistent for the
duration of a project, reusing the same reference for bias correction of all
tumor samples in a cohort.
If your library preparation protocol changes, it's usually best to build a new
reference file and use the new file to analyze the samples prepared under the
new protocol.


Bin-level log2 ratios (.cnr)
----------------------------

In addition to the ``chromosome``, ``start``, ``end``, ``gene``, ``log2`` and
``depth`` columns present in .cnn files, the .cnr file includes each bin's
proportional weight or reliability (``weight``).

The **weight** value is derived from several sources:

- The size of the bin relative to the average bin size (for targets or
  antitargets, separately)
- For a paired or pooled reference, the deviation of the reference log2 value
  from neutral coverage (i.e. distance from 0.0)
- For a pooled reference, the inverse of the variance (i.e. square of ``spread``
  in the reference) of normalized log2 coverage values seen among all normal
  samples at that bin.

This calculated value is used to weight the bin log2 ratio values during
segmentation.
Also, when a genomic region is plotted with CNVkit's "scatter" command, the size
of the plotted datapoints is proportional to each bin's weight -- a relatively
small point indicates a less reliable bin.


Segmented log2 ratios (.cns)
----------------------------

In addition to the ``chromosome``, ``start``, ``end``, ``gene``, ``log2``,
``depth`` and ``weight`` columns present in .cnr files, the .cns file format has
the additional column ``probes``, indicating the number of bins covered by the
segment. CNVkit always writes this column, but treats it as optional on input,
so segment files from other tools can be loaded without it.

The **gene** column concatenates the gene names of all the bins that the segment
covers. The **weight** column sums the bin-level weights, and the **depth** and
**log2** is the weighted mean of the input bin-level values corresponding to
the segment.

Later stages of the pipeline append further columns. The :ref:`segmetrics`
command adds one column per requested statistic, e.g. ``ci_lo`` and ``ci_hi``
for the bounds of the segment mean's confidence interval, and the :ref:`call`
command adds the integer copy number ``cn``.

The **p_ttest** column, added by :ref:`segmetrics` ``--t-test``, is the p-value
of a two-sided one-sample t-test of the segment's constituent bin-level log2
ratios against a mean of 0.0. It thus estimates the probability of observing
bin log2 ratios scattered at least this far from neutral coverage if the
segment were in truth copy-number-neutral, so a small value indicates a segment
whose bins are consistently shifted away from zero. The test is undefined for a
segment covering only one bin, and the field is left empty there. The
``.call.cns`` file written by :ref:`batch` carries this column; the ``.cns``
files written by :ref:`segment` and by ``batch`` do not.
