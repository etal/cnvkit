"""NGS utilities: VCF I/O."""
from __future__ import absolute_import, division, print_function

from ..ngfrills import echo
import collections

def _get_sample(record, sample_id=None):
    """Pick sample by sample ID or defaulting to the first sample.
    """
    if sample_id:
        choices = [s for s in record.samples if s.sample == sample_id]
        assert len(choices) == 1, \
            "Did not find single sample matching sample id %s: %s" \
            % (sample_id, [s.sample for s in record.samples])
        return choices[0]
    else:
        return record.samples[0]


def load_vcf(fname, min_depth=1, sample_id=None,
             skip_hom=True, skip_reject=True, skip_somatic=False):
    """Parse SNV coordinates from a VCF file; group by chromosome.

    Returns a dict of: {chrom: position, zygosity, alt allele frequency}
    """
    chrom_snvs = collections.defaultdict(list)
    for record, zygosity in filter_vcf_lines(fname, min_depth, skip_hom, 
                                             skip_reject, skip_somatic,
                                             sample_id):
        samp = _get_sample(record, sample_id)
        # Read depth
        depth = float(samp.data.DP)
        if "AD" in samp.data._fields:
            # GATK and other callers
            if isinstance(samp.data.AD, (list, tuple)):
                alt_count = float(samp.data.AD[1])
            # VarScan
            else:
                alt_count = float(samp.data.AD)
        elif "AO" in samp.data._fields:
            if samp.data.AO:
                if isinstance(samp.data.AO, (list, tuple)):
                    alt_count = sum([float(x) for x in samp.data.AO])
                else:
                    alt_count = float(samp.data.AO)
            else:
                alt_count = 0
        else:
            echo("Skipping: unsure how to get alternative allele count: %s %s %s %s" %
                 (record.CHROM, record.POS, record.REF, str(samp.data)))
            alt_count = None
        if alt_count is not None:
            alt_freq = alt_count / depth
            chrom_snvs[record.CHROM].append((record.POS, zygosity, alt_freq))
    return chrom_snvs


def filter_vcf_lines(vcf_fname, min_depth, skip_hom, skip_reject, skip_somatic,
                     sample_id=None):
    import vcf
    with open(vcf_fname) as vcffile:
        vcf_reader = vcf.Reader(vcffile)
        for record in vcf_reader:
            if skip_reject and record.FILTER and len(record.FILTER) > 0:
                continue
            if skip_somatic and record.INFO.get("SOMATIC"):
                continue
            # Skip unassigned contigs, alt. HLA haplotypes, etc. XXX
            if len(record.CHROM) > len("chr99"):
                continue
            # Skip homozygous variants (optionally)
            sample = _get_sample(record, sample_id)
            if sample.is_het:
                zygosity = 0.5
            elif sample.gt_type == 0:
                zygosity = 0.0
            else:
                zygosity = 1.0
            if skip_hom and zygosity != 0.5:
                continue
            samp = record.samples[0]
            depth = samp.data.DP
            if depth < min_depth:
                continue
            yield record, zygosity

