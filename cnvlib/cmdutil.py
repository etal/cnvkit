"""Functions reused within command-line implementations."""
from __future__ import absolute_import, division, print_function

import logging

from . import tabio

from .cnary import CopyNumArray as CNA


def read_cna(infile, sample_id=None, meta=None):
    """Read a CNVkit file (.cnn, .cnr, .cns) to create a CopyNumArray object."""
    return tabio.read(infile, into=CNA, sample_id=sample_id, meta=meta)


def load_het_snps(vcf_fname, sample_id, normal_id, min_variant_depth,
                  zygosity_freq, tumor_boost=False):
    if vcf_fname is None:
        return None
    varr = tabio.read(vcf_fname, 'vcf',
                      sample_id=sample_id,
                      normal_id=normal_id,
                      min_depth=min_variant_depth,
                      skip_somatic=True)
    if (zygosity_freq is None and 'n_zygosity' in varr and
        not varr['n_zygosity'].any()):
        # Mutect2 sets all normal genotypes to 0/0 -- work around it
        logging.warn("VCF normal sample's genotypes are all 0/0 or missing; "
                     "inferring genotypes from allele frequency instead")
        zygosity_freq = 0.25
    if zygosity_freq is not None:
        varr = varr.zygosity_from_freq(zygosity_freq, 1 - zygosity_freq)
    if 'n_zygosity' in varr:
        # Infer & drop (more) somatic loci based on genotype
        somatic_idx = (varr['zygosity'] != 0.0) & (varr['n_zygosity'] == 0.0)
        if somatic_idx.any() and not somatic_idx.all():
            logging.info("Skipping %d additional somatic record based on "
                         "T/N genotypes", somatic_idx.sum())
        varr = varr[~somatic_idx]
    orig_len = len(varr)
    varr = varr.heterozygous()
    logging.info("Kept %d heterozygous of %d VCF records",
                 len(varr), orig_len)
    # TODO use/explore tumor_boost option
    if tumor_boost:
        varr['alt_freq'] = varr.tumor_boost()
    return varr


def verify_sample_sex(cnarr, sex_arg, is_male_reference):
    is_sample_female = cnarr.guess_xx(is_male_reference, verbose=False)
    if sex_arg:
        is_sample_female_given = (sex_arg.lower() not in ['y', 'm', 'male'])
        if is_sample_female != is_sample_female_given:
            logging.warn("Sample sex specified as %s "
                         "but chromosomal X/Y ploidy looks like %s",
                         "female" if is_sample_female_given else "male",
                         "female" if is_sample_female else "male")
            is_sample_female = is_sample_female_given
    logging.info("Treating sample %s as %s",
                 cnarr.sample_id or '',
                 "female" if is_sample_female else "male")
    return is_sample_female
