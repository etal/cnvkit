Tumor heterogeneity
===================

DNA samples extracted from solid tumors are rarely completely pure. Stromal or
other normal cells and distinct subclonal tumor-cell populations are typically
present in a sample, and can confound attempts to fit segmented log2 ratio
values to absolute integer copy numbers.

CNVkit provides several points of integration with existing tools and methods
for dealing with tumor heterogeneity and normal-cell contamination.

Estimating tumor purity and normal contamination
------------------------------------------------

A rough estimate of tumor purity can usually be obtained using one or more of
these approaches:

1. A pathologist can visually estimate the purity of an sample taken from a
   solid tumor by examination under a microscope, counting stromal and
   neoplastic cells.
2. If the tumor is belived to be driven by a somatic point mutation, e.g. BRAF
   V600E in melanoma, then that mutation is assumed to be fully clonal and its
   allele frequency indicates the tumor purity. This can be complicated by copy
   number alterations at the same site and whether the point mutation is
   homozygous or heterozygous, but the frequencies of other somatic mutations in
   the same sample may resolve this satisfactorily.
3. Larger-scale, hemizygous losses that cover germline heterozygous SNPs shift
   the allele frequencies of the same SNPs as they are present in the tumor
   sample. In a 50% pure tumor sample, for example, these SNP b-allele
   frequencies would shift from 50% to 67% or 33%, assuming a diploid sample
   (i.e. 1 of 2 copies from the normal sample and 0 or 1 of 1 copy from the
   tumor, depending on whether the variant allele was lost or retained). The
   general calculation is a bit more complicated than in #1 or #2, and can be
   done similarly for copy number gains and homozygous deletions.
4. The log2 ratio values of CNAs in a tumor sample correspond to integer copy
   numbers in tumor cells, and in aggregate these log2 values will cluster
   around values that indicate subclone populations, each with a given ploidy
   and clonality. For example, a single-copy loss in a 50% pure tumor sample
   will have 3/4 the coverage of a neutral site (2/2 normal copies, 1/2 tumor
   copies), for a log2 value of log2(.75) = -0.415. This calculation can also be
   generalized to other copy number states.

Software implementations of the latter three approaches can be used directly on
DNA sequencing data.


Inferring tumor purity and subclonal population fractions from sequencing
-------------------------------------------------------------------------

While inferring the tumor population structure is currently out of the scope of
CNVkit, this work can be done using other third-party programs such as
`THetA2 <http://compbio.cs.brown.edu/projects/theta/>`_,
`PyClone <http://compbio.bccrc.ca/software/pyclone/>`_, or
`BubbleTree <https://www.bioconductor.org/packages/release/bioc/html/BubbleTree.html>`_.
Each of these programs can be used to estimate tumor cell content and infer
integer copy number of tumor subclones in a sample.


Using CNVkit with THetA2
````````````````````````

CNVkit provides wrappers for exporting .cns segments to THetA2's input format
and importing THetA2's result file as CNVkit's segmented .cns files.
See the commands :ref:`export_theta` and :ref:`import-theta` for usage
instructions.

After running the CNVkit :doc:`pipeline` on a sample, and calling SNVs jointly
on the tumor and normal samples, generate the THetA2 input files from the .cns
and .vcf files::

    cnvkit.py export theta Sample_T.cns reference.cnn -v Sample_Paired.vcf

This produces three output files: ``Sample_T.interval_count``,
``Sample_T.tumor.snp_formatted.txt``, and
``Sample_T.normal.snp_formatted.txt``.

Then, run THetA2 (assuming the program was unpacked at ``/path/to/theta2/``)::

    # Generates Sample_T.BEST.results:
    /path/to/theta2/bin/RunTHetA Sample_T.interval_count \
        --TUMOR_SNP Sample_T.tumor.snp_formatted.txt \
        --NORMAL_SNP Sample_T.normal.snp_formatted.txt \
        --BAF --NUM_PROCESSES `nproc` --FORCE

Finally, import THetA2's results back into CNVkit's .cns format, matching the
original segmentation (.cns) to the THetA2-inferred absolute copy number
values.::

    cnvkit.py import-theta Sample_T.cns Sample_T.BEST.results

THetA2 adjusts the segment log2 values to the inferred cellularity of each
detected subclone; this can result in one or two .cns files representing
subclones if more than one clonal tumor cell population was detected. THetA2
also performs some significance testing of each segment representing a CNA, so
there may be fewer segments derived from THetA2 than were originally found by
CNVkit.

The segment values are still log2-transformed in the resulting .cns files, for
convenience in plotting etc. with CNVkit. These files are also easily converted
to other formats using the :ref:`export` command.


Adjusting copy ratios and segments for normal cell contamination
----------------------------------------------------------------

CNVkit's :ref:`call` command uses an estimate of tumor fraction (from
any source) to directly rescale segment log2 ratio values, and SNV b-allele
frequencies if present, to the value that would be seen a completely pure,
uncontaminated sample. Example with tumor purity of 60% and a male reference::

    cnvkit.py call -m none Sample.cns --purity 0.6 -y -o Sample.call.cns

The :ref:`call` command can also convert the segmented log2 ratio estimates to
absolute integer copy numbers. If the tumor cell fraction is known confidently,
use the ``-m clonal`` method to round the log2 ratios to the nearest integer
copy number. Alternatively, the ``-m threshold`` method to applies hard
thresholds. Note that rescaling for purity is optional; either way, integer copy
numbers are emitted unless the ``-m none`` option is used to skip it.

::

    cnvkit.py call -m clonal Sample.cns -y --purity 0.65 -o Sample.call.cns
    # Or, if already rescaled
    cnvkit.py call -m clonal Sample.call.cns -y -o Sample.call.cns
    # With CNVkit's default cutoffs
    cnvkit.py call -m threshold Sample.cns -y -o Sample.call.cns
    # Or, using a custom set of cutoffs
    cnvkit.py call -t=-1.1,-0.4,0.3,0.7 Sample.cns -y -o Sample.call.cns


Export integer copy numbers as BED or VCF
-----------------------------------------

The :ref:`export` ``bed`` and ``vcf`` commands emit integer copy number calls in
the standard BED or VCF formats::

    cnvkit.py export bed Sample.call.cns -y -o Sample.bed
    cnvkit.py export vcf Sample.call.cns -y -o Sample.vcf

If the `.call.cns` files were generated by the :ref:`call` command, the
integer copy numbers calculated in that step will be exported as well.
