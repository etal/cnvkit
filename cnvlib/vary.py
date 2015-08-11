"""NGS utilities: VCF I/O."""
from __future__ import absolute_import, division, print_function

import pandas as pd
import vcf

from .ngfrills import echo


def read_vcf(vcf_fname, sample_id=None, min_depth=1, skip_hom=True,
             skip_reject=False, skip_somatic=True):
    """Parse SNV coordinates from a VCF file into a DataFrame."""
    with open(vcf_fname) as vcffile:
        vcf_reader = vcf.Reader(vcffile)
        rows = parse_records(vcf_reader, sample_id, min_depth, skip_hom,
                             skip_reject, skip_somatic)
        dframe = pd.DataFrame.from_records(rows, columns=[
            "chromosome", "start", "end",
            "ref", "alt", "zygosity", "depth", "alt_count"])
    dframe["alt_freq"] = dframe["alt_count"] / dframe["depth"]
    return dframe


def parse_records(vcf_reader, sample_id, min_depth, skip_hom, skip_reject,
                  skip_somatic):
    """Parse VCF records into DataFrame rows."""
    cnt_reject = 0
    cnt_som = 0
    cnt_depth = 0
    cnt_hom = 0  # DBG
    for record in vcf_reader:
        if skip_reject and record.FILTER and len(record.FILTER) > 0:
            cnt_reject += 1
            continue
        if skip_somatic and record.INFO.get("SOMATIC"):
            cnt_som += 1
            continue

        # Skip unassigned contigs, alt. HLA haplotypes, etc. XXX
        # if len(record.CHROM) > len("chr99"):
        #     continue

        # Skip homozygous variants (optionally)
        sample = _get_sample(record, sample_id)

        # Depth filter
        depth = sample.data.DP
        if depth < min_depth:
            if cnt_depth == 0:
                echo("First rejected depth:", depth)
            cnt_depth += 1
            continue

        # Alt count
        alt_count = _get_alt_count(sample)

        # Het filter
        if sample.is_het:
            zygosity = 0.5
        elif sample.gt_type == 0:
            zygosity = 0.0
        else:
            zygosity = 1.0
        if skip_hom and zygosity != 0.5:
            cnt_hom += 1
            continue

        # Split multiallelics?
        # XXX Ensure sample genotypes are handled properly
        for alt in record.ALT:
            posn = record.POS - 1
            end = _get_end(posn, alt, record.INFO)
            yield (record.CHROM,
                   posn,
                   end,
                   record.REF,
                   str(alt),
                   zygosity,
                   depth,
                   alt_count,
                  )
    echo("Skipped records:", cnt_reject, "reject,",  cnt_som, "somatic,",
         cnt_depth, "depth,", cnt_hom, "homozygous")


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


def _get_alt_count(sample):
    """Get the alternative allele count from a sample in a VCF record."""
    if "AD" in sample.data._fields:
        # GATK and other callers
        if isinstance(sample.data.AD, (list, tuple)):
            alt_count = float(sample.data.AD[1])
        # VarScan
        else:
            alt_count = float(sample.data.AD)
    elif "AO" in sample.data._fields:
        if sample.data.AO:
            if isinstance(sample.data.AO, (list, tuple)):
                alt_count = sum(map(float, sample.data.AO))
            else:
                alt_count = float(sample.data.AO)
        else:
            alt_count = 0
    else:
        echo("Skipping: unsure how to get alternative allele count:",
             sample.site.CHROM, sample.site.POS, sample.site.REF, sample.data)
        alt_count = None  # or 0 or "missing data"?
    return alt_count


def _get_end(posn, alt, info):
    """Get record end position."""
    if "END" in info:
        # Structural variant
        return posn + info['END']
    return posn + len(alt)
