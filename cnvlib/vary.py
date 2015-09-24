"""An array of genomic intervals, treated as variant loci."""
from __future__ import absolute_import, division, print_function

import pandas as pd
import vcf

from . import gary
from .ngfrills import echo


class VariantArray(gary.GenomicArray):
    """An array of genomic intervals, treated as variant loci.

    Required columns: chromosome, start, end, ref, alt
    """
    _required_columns = ("chromosome", "start", "end", "ref", "alt")
    # Extra: zygosity, depth, alt_count, alt_freq

    def __init__(self, data_table, meta_dict=None):
        gary.GenomicArray.__init__(self, data_table, meta_dict)

    # I/O

    @classmethod
    def read_vcf(cls, infile, sample_id=None, min_depth=1, skip_hom=True,
                 skip_reject=False, skip_somatic=True):
        """Parse SNV coordinates from a VCF file into a VariantArray."""
        if isinstance(infile, basestring):
            vcf_reader = vcf.Reader(filename=infile)
        else:
            vcf_reader = vcf.Reader(infile)
        if not vcf_reader.samples:
            raise ValueError("No samples found in VCF file " + str(infile))

        sample_id = _select_sample(vcf_reader, sample_id)
        rows = _parse_records(vcf_reader, sample_id, min_depth, skip_hom,
                             skip_reject, skip_somatic)
        table = pd.DataFrame.from_records(rows, columns=[
            "chromosome", "start", "end", "ref", "alt",
            "zygosity", "depth", "alt_count"])
        table["alt_freq"] = table["alt_count"] / table["depth"]
        return cls(table, {"sample_id": sample_id})

    # def write_vcf(self, outfile=sys.stdout):
    #     """Write data to a file or handle in VCF format."""
    #     # grab from export.export_vcf()


def _select_sample(vcf_reader, sample_id):
    """Select a sample ID in the VCF; ensure it's valid."""
    # ENH - take the paired normal, to select only the tumor records where the
    # normal sample is het
    # def get_mutect_tag(metadata):
    #     for tag in metadata["GATKCommandLine"]:
    #         if tag["ID"] == "MuTect":
    #             return sample_id

    if sample_id is None:
        # Use the VCF header to select the tumor sample
        if "PEDIGREE" in vcf_reader.metadata:
            for tag in vcf_reader.metadata["PEDIGREE"]:
                if "Derived" in tag:
                    sample_id = tag["Derived"]
                    # normal_id = tag["Original"]
                    echo("Selected tumor sample", sample_id,
                         "from the VCF header PEDIGREE tag")
                    break
        elif "GATKCommandLine" in vcf_reader.metadata:
            for tag in vcf_reader.metadata["GATKCommandLine"]:
                if tag.get("ID") == "MuTect":  # any others OK?
                    options = dict(kv.split("=", 1)
                                   for kv in tag["CommandLineOptions"].split()
                                   if '=' in kv)
                    sample_id = options.get('tumor_sample_name')
                    # normal_id = options['normal_sample_name=TR_95_N']
                    echo("Selected tumor sample", sample_id,
                         "from the MuTect VCF header")
                    break

    if sample_id is None:
        sample_id = vcf_reader.samples[0]
    choices = [s for s in vcf_reader.samples if s == sample_id]
    assert len(choices) == 1, (
        "Did not find single sample matching sample id %s: %s"
        % (sample_id, vcf_reader.samples))
    return sample_id


def _parse_records(vcf_reader, sample_id, min_depth, skip_hom, skip_reject,
                  skip_somatic):
    """Parse VCF records into DataFrame rows.

    Apply filters to skip records with low depth, homozygosity, the REJECT
    flag, or the SOMATIC info field.
    """
    cnt_reject = 0  # For logging
    cnt_som = 0
    cnt_depth = 0
    cnt_hom = 0
    for record in vcf_reader:
        if skip_reject and record.FILTER and len(record.FILTER) > 0:
            cnt_reject += 1
            continue
        if skip_somatic and record.INFO.get("SOMATIC"):
            cnt_som += 1
            continue

        sample = record.genotype(sample_id)

        # Depth filter
        if "DP" in sample.data._fields:
            depth = sample.data.DP
        else:
            # SV, probably
            cnt_depth += 1
            continue
        if depth < min_depth:
            cnt_depth += 1
            continue

        # Skip homozygous variants (optionally)
        if sample.is_het:
            zygosity = 0.5
        else:
            if skip_hom:
                cnt_hom += 1
                continue
            if sample.gt_type == 0:
                zygosity = 0.0
            else:
                zygosity = 1.0

        alt_count = _get_alt_count(sample)

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


def _get_alt_count(sample):
    """Get the alternative allele count from a sample in a VCF record."""
    if "AD" in sample.data._fields and sample.data.AD is not None:
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
            alt_count = 0.0
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
