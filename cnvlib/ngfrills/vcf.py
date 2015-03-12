"""NGS utilities: VCF I/O."""
from __future__ import absolute_import, division, print_function

import collections


def load_vcf(fname, min_depth=1, skip_hom=True):
    """Parse SNV coordinates from a VCF file; group by chromosome.

    Returns a dict of: {chrom: position, zygosity, alt allele frequency}
    """
    chrom_snvs = collections.defaultdict(list)
    for record, zygosity in filter_vcf_lines(fname, min_depth, skip_hom):
        # Read depth
        samp = record.samples[0]
        depth = float(samp.data.DP)
        alt_count = float(samp.data.AD[1])
        alt_freq = alt_count / depth
        chrom_snvs[record.CHROM].append((record.POS, zygosity, alt_freq))
    return chrom_snvs


def filter_vcf_lines(vcf_fname, min_depth, skip_hom):
    import vcf
    with open(vcf_fname) as vcffile:
        vcf_reader = vcf.Reader(vcffile)
        for record in vcf_reader:
            # Skip unassigned contigs
            if len(record.CHROM) > len("chr22"):
                continue
            # Skip homozygous variants (optionally)
            # XXX if no 'AF', check 'FA' (for MuTect)
            zygosity = record.INFO['AF'][0]  # 1.0 or 0.5
            if skip_hom and zygosity != 0.5:
                continue
            samp = record.samples[0]
            depth = samp.data.DP
            if depth < min_depth:
                continue
            yield record, zygosity

