"""Functions reused within command-line implementations."""
from __future__ import absolute_import, division, print_function

import logging
import sys

from skgenome import tabio

from .cnary import CopyNumArray as CNA


def read_cna(infile, sample_id=None, meta=None):
    """Read a CNVkit file (.cnn, .cnr, .cns) to create a CopyNumArray object."""
    return tabio.read(infile, into=CNA, sample_id=sample_id, meta=meta)


def load_het_snps(vcf_fname, sample_id=None, normal_id=None,
                  min_variant_depth=20, zygosity_freq=None, tumor_boost=False):
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
        logging.warning("VCF normal sample's genotypes are all 0/0 or missing; "
                        "inferring genotypes from allele frequency instead")
        zygosity_freq = 0.25
    if zygosity_freq is not None:
        varr = varr.zygosity_from_freq(zygosity_freq, 1 - zygosity_freq)
    if 'n_zygosity' in varr:
        # Infer & drop (more) somatic loci based on genotype
        somatic_idx = (varr['zygosity'] != 0.0) & (varr['n_zygosity'] == 0.0)
        if somatic_idx.any() and not somatic_idx.all():
            logging.info("Skipping %d additional somatic records based on "
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
            logging.warning("Sample sex specified as %s "
                            "but chromosomal X/Y ploidy looks like %s",
                            "female" if is_sample_female_given else "male",
                            "female" if is_sample_female else "male")
            is_sample_female = is_sample_female_given
    logging.info("Treating sample %s as %s",
                 cnarr.sample_id or '',
                 "female" if is_sample_female else "male")
    return is_sample_female


def write_tsv(outfname, rows, colnames=None):
    """Write rows, with optional column header, to tabular file."""
    with tabio.safe_write(outfname or sys.stdout) as handle:
        if colnames:
            header = '\t'.join(colnames) + '\n'
            handle.write(header)
        handle.writelines('\t'.join(map(str, row)) + '\n'
                           for row in rows)


def write_text(outfname, text, *more_texts):
    """Write one or more strings (blocks of text) to a file."""
    with tabio.safe_write(outfname or sys.stdout) as handle:
        handle.write(text)
        if more_texts:
            for mtext in more_texts:
                handle.write(mtext)


def write_dataframe(outfname, dframe, header=True):
    """Write a pandas.DataFrame to a tabular file."""
    with tabio.safe_write(outfname or sys.stdout) as handle:
        dframe.to_csv(handle, header=header,
                      index=False, sep='\t', float_format='%.6g')
