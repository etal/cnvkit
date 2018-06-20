#!/bin/bash
# Call the given germline-heterozygous SNPs in a tumor sample
# Output filename is "calls.vcf"
# NB: Chromosome names must match between the BAM, VCF, and (optional) FASTA.
# Without the FASTA genome sequence, the VCF is a little weird but still usable.
# If given, FASTA should be bgzipped and indexed with 'samtools faidx'.

if [ -z "$1" ]
then
    echo "usage: $0 bam vcf [fasta]"
    exit 1
fi

bam=$1
vcf=$2
fasta=$3
# TODO make sure chroms match

if [ -z "$fasta" ]
then
    fa_opts=""
else
    fa_opts="--fasta-ref $fasta"
    bgzip -r $fasta
    if [ ! -e "${fasta}.fai" ]
    then
        samtools faidx $fasta
    fi
fi

# TODO/ENH: Can bcftools emit just the sites, no samples? Should it?
bcftools view --exclude-uncalled --trim-alt-alleles --genotype ^miss \
    --output-type v --output-file _snps.vcf $vcf
grep -v '^#' _snps.vcf | cut -f1,2 > sites.txt
samtools mpileup \
    --ignore-RG --skip-indels --count-orphans --output-tags DP,AD \
    $fa_opts \
    --positions sites.txt \
    --VCF --uncompressed --output calls.vcf \
    $bam
