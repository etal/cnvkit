"""An array of genomic intervals, treated as variant loci."""
from __future__ import absolute_import, division, print_function
from builtins import next
from builtins import str
from builtins import map
from builtins import zip
from past.builtins import basestring

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
    # Extra: somatic, zygosity, depth, alt_count, alt_freq

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
    def read_vcf(cls, infile, sample_id=None, normal_id=None, min_depth=None,
                 skip_hom=False, skip_reject=False, skip_somatic=False):
        """Parse SNV coordinates from a VCF file into a VariantArray."""
        if isinstance(infile, basestring):
            vcf_reader = vcf.Reader(filename=infile)
        else:
            vcf_reader = vcf.Reader(infile)
        if not vcf_reader.samples:
            logging.warn("VCF file %s has no samples; parsing minimal info",
                         infile)
            return cls._read_vcf_nosample(infile, sample_id, skip_reject)

        columns = [
            "chromosome", "start", "end", "ref", "alt",
            "somatic", "zygosity", "depth", "alt_count"]
        sample_id, normal_id = _select_sample(vcf_reader, sample_id, normal_id)
        if normal_id:
            columns.extend(["n_zygosity", "n_depth", "n_alt_count"])
        rows = _parse_records(vcf_reader, sample_id, normal_id, skip_reject)
        table = pd.DataFrame.from_records(rows, columns=columns)
        table["alt_freq"] = table["alt_count"] / table["depth"]
        if normal_id:
            table["n_alt_freq"] = table["n_alt_count"] / table["n_depth"]
        # Filter out records as requested
        cnt_depth = cnt_hom = cnt_som = 0
        if min_depth:
            dkey = "n_depth" if "n_depth" in table else "depth"
            idx_depth = table[dkey] >= min_depth
            cnt_depth = (~idx_depth).sum()
            table = table[idx_depth]
        if skip_hom:
            zkey = "n_zygosity" if "n_zygosity" in table else "zygosity"
            idx_het = (table[zkey] != 0.0) & (table[zkey] != 1.0)
            cnt_hom = (~idx_het).sum()
            table = table[idx_het]
        if skip_somatic:
            idx_som = table["somatic"]
            cnt_som = idx_som.sum()
            table = table[~idx_som]
        logging.info("Skipped records: %d somatic, %d depth, %d homozygous",
                     cnt_som, cnt_depth, cnt_hom)

        return cls(table, {"sample_id": sample_id})

    @classmethod
    def _read_vcf_nosample(cls, vcf_file, sample_id=None, skip_reject=False):
        table = pd.read_table(vcf_file,
                              comment="#",
                              header=None,
                              na_filter=False,
                              names=["chromosome", "start", "_ID", "ref", "alt",
                                     "_QUAL", "filter", "info"],
                              usecols=cls._required_columns,
                              # usecols=["chromosome", "start", "ref", "alt",
                              #          # "filter", "info",
                              #         ],
                              # ENH: converters=func -> to parse each col
                              dtype=dict(list(zip(cls._required_columns,
                                             cls._required_dtypes))),
                             )
        # ENH: do things with filter, info
        # if skip_reject and record.FILTER and len(record.FILTER) > 0:
        table['end'] = table['start']  # TODO: _get_end
        table['start'] -= 1
        table = table.loc[:, cls._required_columns]
        return cls(table, {"sample_id": sample_id})


    # def write_vcf(self, outfile=sys.stdout):
    #     """Write data to a file or handle in VCF format."""
    #     # grab from export.export_vcf()


def _select_sample(vcf_reader, sample_id, normal_id):
    """Select a sample ID in the VCF; ensure it's valid."""
    peds = list(_parse_pedigrees(vcf_reader))
    if sample_id is None and normal_id is None and peds:
        # Use the VCF header to select the tumor sample
        # Take the first entry, if any
        sample_id, normal_id = peds[0]
    elif sample_id and normal_id:
        # Take the given IDs at face value, just verify below
        pass
    elif sample_id:
        if peds:
            for sid, nid in peds:
                if sid == sample_id:
                    normal_id = nid
                    break
        # if normal_id is None and len(vcf_reader.samples == 2):
        #     normal_id = next(s for s in vcf_reader.samples if s != sample_id)
    elif normal_id:
        if peds:
            for sid, nid in peds:
                if nid == normal_id:
                    sample_id = sid
                    break
        else:
            try:
                sample_id = next(s for s in vcf_reader.samples if s != normal_id)
            except StopIteration:
                raise ValueError(
                    "No other sample in VCF besides the specified normal " +
                    normal_id + "; did you mean to use this as the sample_id "
                    "instead?")
    else:
        sample_id = vcf_reader.samples[0]

    _confirm_unique(sample_id, vcf_reader.samples)
    if normal_id:
        _confirm_unique(normal_id, vcf_reader.samples)
    logging.info("Selected test sample " + sample_id +
                 (" and control sample %s" % normal_id if normal_id else ''))
    return sample_id, normal_id


def _parse_pedigrees(vcf_reader):
    """Extract tumor/normal pair sample IDs from the VCF header.

    Return an iterable of (tumor sample ID, normal sample ID).
    """
    if "PEDIGREE" in vcf_reader.metadata:
        for tag in vcf_reader.metadata["PEDIGREE"]:
            if "Derived" in tag:
                sample_id = tag["Derived"]
                normal_id = tag["Original"]
                logging.debug("Found tumor sample %s and normal sample %s "
                              "in the VCF header PEDIGREE tag",
                              sample_id, normal_id)
                yield sample_id, normal_id

    elif "GATKCommandLine" in vcf_reader.metadata:
        for tag in vcf_reader.metadata["GATKCommandLine"]:
            if tag.get("ID") == "MuTect":  # any others OK?
                options = dict(kv.split("=", 1)
                                for kv in tag["CommandLineOptions"].split()
                                if '=' in kv)
                sample_id = options.get('tumor_sample_name')
                normal_id = options['normal_sample_name']
                logging.debug("Found tumor sample %s and normal sample "
                              "%s in the MuTect VCF header",
                              sample_id, normal_id)
                yield sample_id, normal_id


def _confirm_unique(sample_id, samples):
    occurrences = [s for s in samples if s == sample_id]
    if len(occurrences) != 1:
        raise ValueError(
            "Did not find a single sample ID '%s' in: %s"
            % (sample_id, samples))


def _parse_records(vcf_reader, sample_id, normal_id, skip_reject):
    """Parse VCF records into DataFrame rows.

    Apply filters to skip records with low depth, homozygosity, the REJECT
    flag, or the SOMATIC info field.
    """
    cnt_reject = 0  # For logging
    for record in vcf_reader:
        is_som = False
        if skip_reject and record.FILTER and len(record.FILTER) > 0:
            cnt_reject += 1
            continue
        if record.INFO.get("SOMATIC"):
            is_som = True

        sample = record.genotype(sample_id)
        depth, zygosity, alt_count = _extract_genotype(sample)
        if normal_id:
            normal = record.genotype(normal_id)
            n_depth, n_zygosity, n_alt_count = _extract_genotype(normal)
            if n_zygosity == 0:
                is_som = True

        # Split multiallelics?
        # XXX Ensure sample genotypes are handled properly
        for alt in record.ALT:
            posn = record.POS - 1
            end = _get_end(posn, alt, record.INFO)
            row = (record.CHROM, posn, end, record.REF, str(alt),
                   is_som,
                   zygosity,
                   depth,
                   alt_count,
                  )
            if normal_id:
                row += (n_zygosity, n_depth, n_alt_count)
            yield row

    if cnt_reject:
        logging.info('Filtered out %d records', cnt_reject)


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
        return info['END']
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
    seg_bafs = np.asarray([(np.median(subvarr.mirrored_baf())
                            if len(subvarr) else np.nan)
                           for _seg, subvarr in varr.by_ranges(segarr)])
    cn1 = 0.5 * (1 - seg_bafs) * seg_depths
    cn2 = seg_depths - cn1
    # segout = segarr.copy()
    # segout.update({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
    # return segout
    return pd.DataFrame({"baf": seg_bafs, "CN1": cn1, "CN2": cn2})
