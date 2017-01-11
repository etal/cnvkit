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

Estimation from a SNP b-allele frequencies works by comparing the shift in
allele frequency of heterozygous, germline SNPs in the tumor sample from the
expected ~50% -- e.g. a 3-copy segment in a diploid genome would have log2
ratio of +0.58 and heterozygous SNPs would have an average BAF of 67% or 33% if
the tumor sample is fully clonal, and closer to log2 0.0 and BAF 0.5 if there
is normal-cell contamination and/or :doc:`tumor heterogeneity <heterogeneity>`.

Typically you would use a properly formatted :ref:`vcfformat` from joint
tumor-normal SNV calling, e.g. the output of MuTect, VarDict, or FreeBayes,
having already flagged somatic mutations so they can be skipped in this
analysis. If you have no matched normal sample for a given tumor, you can use
1000 Genomes common SNP sites to extract the likely germline SNVs from a
tumor-only VCF, and use just those sites with THetA2 (or another tool like
PyClone or BubbleTree).

Use SNP b-allele frequencies from a VCF in these commands:

- :ref:`call`
- :ref:`scatter`
- :ref:`export` ``nexus-ogt`` and ``theta``

