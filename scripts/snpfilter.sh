#!/bin/bash
# Call the given germline-heterozygous SNPs in a tumor sample
# Output filename is "calls.vcf"
# NB: Chromosome names must match between the BAM, VCF, and (optional) FASTA.
# Without the FASTA genome sequence, the VCF is a little weird but still usable.
# If given, FASTA should be bgzipped and indexed with 'samtools faidx'.

bam=$1
vcf=$2
# TODO make fasta optional & handle safely
#fasta=$3
# TODO make sure chroms match

# Can bcftools emit just the sites, no samples? should it? will mpileup accept?
bcftools view --exclude-uncalled --trim-alt-alleles --genotype ^miss \
    --output-type v --output-file _snps.vcf $vcf
grep -v '^#' _snps.vcf | cut -f1,2 > sites.txt
bgzip -r $fasta
samtools mpileup \
    --ignore-RG --skip-indels --count-orphans \
    --output-tags DP,AD \
    --fasta-ref $fasta \
    --positions sites.txt \
    --VCF --uncompressed --output calls.vcf \
    $bam

