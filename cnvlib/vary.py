"""An array of genomic intervals, treated as variant loci."""
from __future__ import absolute_import, division, print_function

import logging

import numpy as np
import pandas as pd
import vcf

from . import gary


class VariantArray(gary.GenomicArray):
    """An array of genomic intervals, treated as variant loci.

    Required columns: chromosome, start, end, ref, alt
    """
    _required_columns = ("chromosome", "start", "end", "ref", "alt")
    _required_dtypes = ("string", "int", "int", "string", "string")
    # Extra: zygosity, depth, alt_count, alt_freq

    def __init__(self, data_table, meta_dict=None):
        gary.GenomicArray.__init__(self, data_table, meta_dict)

    def mirrored_baf(self, above_half=True, tumor_boost=True):
        """Mirrored B-allele frequencies (BAFs).

        Flip BAFs to be all above 0.5 (if `above_half`) or below 0.5, for
        consistency.

        With `tumor_boost`, normalize tumor-sample allele frequencies to the
        matched normal allele frequencies.
        """
        if tumor_boost and "n_alt_freq" in self:
            alt_freqs = self.tumor_boost()
        else:
            alt_freqs = self["alt_freq"]
        fromhalf = np.abs(alt_freqs - 0.5)
        if above_half:
            return fromhalf + 0.5
        else:
            return 0.5 - fromhalf

    def tumor_boost(self):
        """TumorBoost normalization of tumor-sample allele frequencies.

        De-noises the signal for detecting LOH.
        """
        if "n_alt_freq" in self:
            n_freqs = self["n_alt_freq"]
        else:
            raise ValueError("TumorBoost requires a paired normal sample in "
                             "the VCF.")
        t_freqs = self["alt_freq"]
        return _tumor_boost(t_freqs, n_freqs)

    # I/O

    @classmethod
    def read_vcf(cls, infile, sample_id=None, normal_id=None, min_depth=1,
                 skip_hom=True, skip_reject=False, skip_somatic=True):
        """Parse SNV coordinates from a VCF file into a VariantArray."""
        if isinstance(infile, basestring):
            vcf_reader = vcf.Reader(filename=infile)
        else:
            vcf_reader = vcf.Reader(infile)
        if not vcf_reader.samples:
            raise ValueError("No samples found in VCF file " + str(infile))

        sample_id, normal_id = _select_sample(vcf_reader, sample_id, normal_id)
        if normal_id:
            columns = [
                "chromosome", "start", "end", "ref", "alt",
                "zygosity", "depth", "alt_count",
                "n_zygosity", "n_depth", "n_alt_count"]
        else:
            columns = [
                "chromosome", "start", "end", "ref", "alt",
                "zygosity", "depth", "alt_count"]
        rows = _parse_records(vcf_reader, sample_id, normal_id, min_depth,
                              skip_hom, skip_reject, skip_somatic)
        table = pd.DataFrame.from_records(rows, columns=columns)
        table["alt_freq"] = table["alt_count"] / table["depth"]
        if normal_id:
            table["n_alt_freq"] = table["n_alt_count"] / table["n_depth"]
        return cls(table, {"sample_id": sample_id})


    # def write_vcf(self, outfile=sys.stdout):
    #     """Write data to a file or handle in VCF format."""
    #     # grab from export.export_vcf()


def _select_sample(vcf_reader, sample_id, normal_id):
    """Select a sample ID in the VCF; ensure it's valid."""
    if sample_id is None and normal_id is None:
        # Use the VCF header to select the tumor sample
        if "PEDIGREE" in vcf_reader.metadata:
            for tag in vcf_reader.metadata["PEDIGREE"]:
                if "Derived" in tag:
                    sample_id = tag["Derived"]
                    normal_id = tag["Original"]
                    logging.info("Selected tumor sample %s and normal sample %s"
                                 " from the VCF header PEDIGREE tag",
                                 sample_id, normal_id)
                    break
        elif "GATKCommandLine" in vcf_reader.metadata:
            for tag in vcf_reader.metadata["GATKCommandLine"]:
                if tag.get("ID") == "MuTect":  # any others OK?
                    options = dict(kv.split("=", 1)
                                   for kv in tag["CommandLineOptions"].split()
                                   if '=' in kv)
                    sample_id = options.get('tumor_sample_name')
                    normal_id = options['normal_sample_name']
                    logging.info("Selected tumor sample %s and normal sample "
                                 "%s from the MuTect VCF header",
                                 sample_id, normal_id)
                    break

    if sample_id:
        pass
        # if normal_id is None and len(vcf_reader.samples == 2):
        #     normal_id = next(s for s in vcf_reader.samples if s != sample_id)
    elif normal_id:
        try:
            sample_id = next(s for s in vcf_reader.samples if s != normal_id)
        except StopIteration:
            raise ValueError("No other sample in VCF besides the specified " +
                             "normal " + normal_id + "; did you mean to use" +
                             "this as the sample_id instead?")
    else:
        sample_id = vcf_reader.samples[0]

    _confirm_unique(sample_id, vcf_reader.samples)
    if normal_id:
        _confirm_unique(normal_id, vcf_reader.samples)
    return sample_id, normal_id


def _confirm_unique(sample_id, samples):
    occurrences = [s for s in samples if s == sample_id]
    if len(occurrences) != 1:
        raise ValueError(
            "Did not find a single sample ID '%s' in: %s"
            % (sample_id, samples))


def _parse_records(vcf_reader, sample_id, normal_id, min_depth,
                   skip_hom, skip_reject, skip_somatic):
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
        depth, zygosity, alt_count = _extract_genotype(sample)
        if normal_id:
            normal = record.genotype(normal_id)
            n_depth, n_zygosity, n_alt_count = _extract_genotype(normal)
            if n_depth is None or n_depth < min_depth:
                cnt_depth += 1
                continue
            if skip_hom and n_zygosity in (0.0, 1.0):
                cnt_hom += 1
                continue
            if skip_somatic and n_zygosity == 0:
                cnt_som += 1
                continue
        else:
            if depth is None or depth < min_depth:
                cnt_depth += 1
                continue
            if skip_hom and zygosity  in (0.0, 1.0):
                cnt_hom += 1
                continue

        # Split multiallelics?
        # XXX Ensure sample genotypes are handled properly
        for alt in record.ALT:
            posn = record.POS - 1
            end = _get_end(posn, alt, record.INFO)
            row = (record.CHROM, posn, end, record.REF, str(alt),
                zygosity,
                depth,
                alt_count,
                )
            if normal_id:
                row += (n_zygosity, n_depth, n_alt_count)
            yield row

    logging.info("Skipped records: %d reject, %d somatic, %d depth, "
                 "%d homozygous", cnt_reject, cnt_som, cnt_depth, cnt_hom)


def _extract_genotype(sample):
    if "DP" in sample.data._fields:
        depth = sample.data.DP
    else:
        # SV, probably
        depth = alt_count = None
    if sample.is_het:
        zygosity = 0.5
    elif sample.gt_type == 0:
        zygosity = 0.0
    else:
        zygosity = 1.0
    alt_count = (_get_alt_count(sample) if sample.gt_type else 0.0)
    return depth, zygosity, alt_count


def _get_alt_count(sample):
    """Get the alternative allele count from a sample in a VCF record."""
    if "AD" in sample.data._fields and sample.data.AD is not None:
        # GATK and other callers
        if isinstance(sample.data.AD, (list, tuple)):
            alt_count = float(sample.data.AD[1])
        # VarScan
        else:
            alt_count = float(sample.data.AD)
    elif "CLCAD2" in sample.data._fields and sample.data.CLCAD2 is not None:
        # Qiagen CLC Genomics Server -- similar to GATK's AD
        alt_count = float(sample.data.CLCAD2[1])
    elif "AO" in sample.data._fields:
        if sample.data.AO:
            if isinstance(sample.data.AO, (list, tuple)):
                alt_count = sum(map(float, sample.data.AO))
            else:
                alt_count = float(sample.data.AO)
        else:
            alt_count = 0.0
    else:
        logging.warn("Skipping %s:%d %s; "
                     "unsure how to get alternative allele count: %s",
                     sample.site.CHROM, sample.site.POS, sample.site.REF,
                     sample.data)
        alt_count = None  # or 0 or "missing data"?
    return alt_count


def _get_end(posn, alt, info):
    """Get record end position."""
    if "END" in info:
        # Structural variant
        return posn + info['END']
    return posn + len(alt)


def _tumor_boost(t_freqs, n_freqs):
    """Normalize tumor-sample allele frequencies.

    boosted = { 0.5 (t/n)           if t < n
                1 - 0.5(1-t)/(1-n)  otherwise

    See: TumorBoost, Bengtsson et al. 2010
    """
    lt_mask = (t_freqs < n_freqs)
    lt_idx = np.nonzero(lt_mask)[0]
    gt_idx = np.nonzero(~lt_mask)[0]
    out = np.zeros_like(t_freqs)
    out[lt_idx] = 0.5 * t_freqs.take(lt_idx) / n_freqs.take(lt_idx)
    out[gt_idx] = 1 - 0.5 * (1 - t_freqs.take(gt_idx)
                            ) / (1 - n_freqs.take(gt_idx))
    return out


def _allele_specific_copy_numbers(segarr, varr, ploidy=2):
    """Split total copy number between alleles based on BAF.

    See: PSCBS, Bentsson et al. 2011
    """
    seg_depths = ploidy * np.exp2(segarr["log2"])
    seg_bafs = [np.median(subvarr.mirrored_baf())
                for _seg, subvarr in varr.by_ranges(segarr)]
    cn1 = 0.5 * (1 - seg_bafs) * seg_depths
    cn2 = seg_depths - cn1
    # segout = segarr.copy()
    # segout.update({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
    # return segout
    return pd.DataFrame({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
