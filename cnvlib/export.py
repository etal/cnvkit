"""Export CNVkit objects and files to other formats."""
from __future__ import absolute_import, division, print_function

import collections
import logging

import numpy as np
import pandas as pd
from Bio._py3k import map, range, zip

from . import call, core
from .cnary import CopyNumArray as CNA
from .vary import VariantArray as VA


def merge_samples(filenames):
    """Merge probe values from multiple samples into a 2D table (of sorts).

    Input:
        dict of {sample ID: (probes, values)}
    Output:
        list-of-tuples: (probe, log2 coverages...)
    """
    def label_with_gene(cnarr):
        row2label = lambda row: "{}:{}-{}:{}".format(
            row['chromosome'], row['start'], row['end'], row['gene'])
        return cnarr.data.apply(row2label, axis=1)

    if not filenames:
        return []
    first_cnarr = CNA.read(filenames[0])
    out_table = first_cnarr.data.loc[:, ["chromosome", "start", "end", "gene"]]
    out_table["label"] = label_with_gene(first_cnarr)
    out_table[first_cnarr.sample_id] = first_cnarr["log2"]
    for fname in filenames[1:]:
        cnarr = CNA.read(fname)
        # Verify labels match
        labels = label_with_gene(cnarr)
        if not (labels == out_table["label"]).all():
            raise ValueError("Mismatched row coordinates in %s" % fname)
        # Copy the next column by sample ID
        if cnarr.sample_id in out_table.columns:
            raise ValueError("Duplicate sample ID: %s" % cnarr.sample_id)
        out_table[cnarr.sample_id] = cnarr["log2"]
        del cnarr
    return out_table


# Supported formats:

def fmt_cdt(sample_ids, table):
    """Format as CDT."""
    outheader = ['GID', 'CLID', 'NAME', 'GWEIGHT'] + sample_ids
    header2 = ['AID', '', '', '']
    header2.extend(['ARRY' + str(i).zfill(3) + 'X'
                    for i in range(len(sample_ids))])
    outrows = [header2]
    outtable = pd.concat([
        pd.DataFrame({
            "GID": table.index.apply(lambda x: "GENE%dX" % x),
            "CLID": table.index.apply(lambda x: "IMAGE:%d" % x),
            "NAME": table["label"],
            "GWEIGHT": 1,
        }),
        table.drop(["chromosome", "start", "end", "gene", "label"],
                   axis=1)],
        axis=1)
    outrows.extend(outtable.itertuples(index=False))
    return outheader, outrows


# TODO
def fmt_gct(sample_ids, table):
    return NotImplemented


def fmt_jtv(sample_ids, table):
    """Format for Java TreeView."""
    outheader = ["CloneID", "Name"] + sample_ids
    outtable = pd.concat([
        pd.DataFrame({
            "CloneID": "IMAGE:",
            "Name": table["label"],
        }),
        table.drop(["chromosome", "start", "end", "gene", "label"],
                   axis=1)],
        axis=1)
    outrows = outtable.itertuples(index=False)
    return outheader, outrows


# Special cases

def export_nexus_basic(sample_fname):
    """Biodiscovery Nexus Copy Number "basic" format.

    Only represents one sample per file.
    """
    cnarr = CNA.read(sample_fname)
    out_table = cnarr.data.loc[:, ['chromosome', 'start', 'end', 'gene', 'log2']]
    out_table['probe'] = cnarr.labels()
    return out_table


def export_nexus_ogt(sample_fname, vcf_fname):
    """Biodiscovery Nexus Copy Number "Custom-OGT" format.

    To create the b-allele frequencies column, alterate allele frequencies from
    the VCF are aligned to the .cnr file bins.  Bins that contain no variants
    are left blank; if a bin contains multiple variants, then the frequencies
    are all "mirrored" to be above .5, then the median of those values is taken.
    """
    def mirrored_baf_median(vals):
        shift = np.median(np.abs(vals - .5))
        if np.median(vals) > .5:
            return .5 + shift
        else:
            return .5 - shift

    cnarr = CNA.read(sample_fname)
    varr = VA.read_vcf(vcf_fname)
    bafs = cnarr.match_to_bins(varr, 'alt_freq', np.nan,
                               summary_func=mirrored_baf_median)
    logging.info("Placed %d variants into %d bins",
                 sum(~np.isnan(bafs)), len(cnarr))
    out_table = cnarr.data.loc[:, ['chromosome', 'start', 'end', 'log2']]
    out_table = out_table.rename(columns={
        "chromosome": "Chromosome",
        "start": "Position",
        "end": "Position",
        "log2": "Log R Ratio",
    })
    out_table["B-Allele Frequency"] = bafs
    return out_table


def export_seg(sample_fnames):
    """SEG format for copy number segments.

    Segment breakpoints are not the same across samples, so samples are listed
    in serial with the sample ID as the left column.
    """
    out_tables = []
    chrom_ids = None
    for fname in sample_fnames:
        segments = CNA.read(fname)
        if chrom_ids is None:
            # Create & store
            chrom_ids = create_chrom_ids(segments)
        else:
            # Verify
            core.assert_equal("Segment chromosome names differ",
                              previous=chrom_ids.keys(),
                              current=create_chrom_ids(segments).keys())
        table = segments.data.loc[:, ["start", "end"]]
        table["ID"] = segments.sample_id
        table["mean"] = segments.data["log2"]
        table["chromosome"] = [chrom_ids[chrom]
                               for chrom in segments["chromosome"]]
        if "probes" in segments:
            table["num_probes"] = segments["probes"]
            sorted_cols = ["ID", "chromosome", "start", "end", "num_probes",
                           "mean"]
        else:
            sorted_cols = ["ID", "chromosome", "start", "end", "mean"]
        out_tables.append(table.reindex(columns=sorted_cols))
    return pd.concat(out_tables)


def create_chrom_ids(segments):
    """Map chromosome names to integers in the order encountered."""
    mapping = collections.OrderedDict()
    curr_idx = 1
    for chrom in segments.chromosome:
        if chrom not in mapping:
            mapping[chrom] = curr_idx
            curr_idx += 1
    return mapping


# _____________________________________________________________________________
# BED

def export_bed(sample_fnames, ploidy, is_reference_male,
               sample_id=None, show="ploidy"):
    """Export to BED format.

    For each region in each sample (possibly filtered according to `show`),
    the columns are:

        - reference sequence name
        - start (0-indexed)
        - end
        - sample name or given label
        - integer copy number

    By default (show="ploidy"), skip regions where copy number is the default
    ploidy, i.e. equal to 2 or the value set by --ploidy.
    If show="variant", skip regions where copy number is neutral, i.e. equal to
    the reference ploidy on autosomes, or half that on sex chromosomes.
    """
    bed_tables = []
    for fname in sample_fnames:
        segs = CNA.read(fname)
        tbl = segments2bed(segs, sample_id or segs.sample_id, ploidy,
                           is_reference_male, show)
        bed_tables.append(tbl)
    return pd.concat(bed_tables)


def segments2bed(segments, label, ploidy, is_reference_male, show):
    """Convert a copy number array to a BED-like DataFrame."""
    absolutes = call.absolute_pure(segments, ploidy, is_reference_male)
    out = segments.data.loc[:, ["chromosome", "start", "end"]]
    out["label"] = label
    out["ncopies"] = np.rint(absolutes)
    if show == "ploidy":
        # Skip regions of default ploidy
        out = out[out["ncopies"] != ploidy]
    elif show == "variants":
        # Skip regions of non-neutral copy number
        # XXX: female samples w/ male reference will show chrX as variant
        haploid_mask = (segments["chromsome"] == segments._chr_y_label)
        if is_reference_male:
            haploid_mask |= (segments["chromsome"] == segments._chr_x_label)
        out = out[(~haploid_mask & (out["ncopies"] != ploidy)) |
                  (haploid_mask & (out["ncopies"] != ploidy//2))]
    return out


# _____________________________________________________________________________
# VCF

VCF_HEADER = """\
##fileformat=VCFv4.0
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
"""
# #CHROM  POS   ID  REF ALT   QUAL  FILTER  INFO  FORMAT  NA00001
# 1 2827693   . CCGTGGATGCGGGGACCCGCATCCCCTCTCCCTTCACAGCTGAGTGACCCACATCCCCTCTCCCCTCGCA  C . PASS  SVTYPE=DEL;END=2827680;BKPTID=Pindel_LCS_D1099159;HOMLEN=1;HOMSEQ=C;SVLEN=-66 GT:GQ 1/1:13.9
# 2 321682    . T <DEL>   6 PASS    IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-105;CIPOS=-56,20;CIEND=-10,62  GT:GQ 0/1:12
# 3 12665100  . A <DUP>   14  PASS  IMPRECISE;SVTYPE=DUP;END=12686200;SVLEN=21100;CIPOS=-500,500;CIEND=-500,500   GT:GQ:CN:CNQ  ./.:0:3:16.2
# 4 18665128  . T <DUP:TANDEM>  11  PASS  IMPRECISE;SVTYPE=DUP;END=18665204;SVLEN=76;CIPOS=-10,10;CIEND=-10,10  GT:GQ:CN:CNQ  ./.:0:5:8.3


def export_vcf(segments, ploidy, is_reference_male, is_sample_female,
               sample_id=None):
    """Convert segments to Variant Call Format.

    For now, only 1 sample per VCF. (Overlapping CNVs seem tricky.)

    Spec: https://samtools.github.io/hts-specs/VCFv4.2.pdf
    """
    vcf_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                   "INFO", "FORMAT", sample_id or segments.sample_id]
    vcf_rows = segments2vcf(segments, ploidy, is_reference_male,
                            is_sample_female)
    table = pd.DataFrame.from_records(vcf_rows, columns=vcf_columns)
    vcf_body = table.to_csv(sep='\t', header=True, index=False,
                            float_format="%.3g")
    return VCF_HEADER, vcf_body


def segments2vcf(segments, ploidy, is_reference_male, is_sample_female):
    """Convert copy number segments to VCF records."""
    out_dframe = segments.data.loc[:, ["chromosome", "end", "log2", "probes"]]
    abs_dframe = call.absolute_dataframe(segments, ploidy, 1.0,
                                         is_reference_male, is_sample_female)
    out_dframe["ncopies"] = np.rint(abs_dframe["absolute"])
    idx_losses = (out_dframe["ncopies"] < abs_dframe["expect"])

    starts = segments.start.copy()
    starts[starts == 0] = 1
    out_dframe["start"] = starts

    svlen = segments.end - segments.start
    svlen[idx_losses] *= -1
    out_dframe["svlen"] = svlen

    out_dframe["svtype"] = "DUP"
    out_dframe.loc[idx_losses, "svtype"] = "DEL"

    out_dframe["format"] = "GT:GQ:CN:CNQ"
    out_dframe.loc[idx_losses, "format"] = "GT:GQ" # :CN:CNQ ?

    # Reformat this data to create INFO and genotype
    # TODO be more clever about this
    for (_idx, out_row), (_idx, abs_row) in zip(out_dframe.iterrows(),
                                                abs_dframe.iterrows()):
        if out_row["ncopies"] == abs_row["expect"] or not str(out_row["probes"]).isdigit():
            # Skip regions of neutral copy number
            continue  # or "CNV" for subclonal?

        if out_row["ncopies"] > abs_row["expect"]:
            genotype = "0/1:0:%d:%g" % (out_row["ncopies"], int(out_row["probes"]))
        elif out_row["ncopies"] < abs_row["expect"]:
            # TODO XXX handle non-diploid ploidies, haploid chroms
            if out_row["ncopies"] == 0:
                # Complete deletion, 0 copies
                gt = "1/1"
            else:
                # Single copy deletion
                gt = "0/1"
            genotype = "%s:%d" % (gt, int(out_row["probes"]))

        info = ";".join(["IMPRECISE",
                         "SVTYPE=%s" % out_row["svtype"],
                         "END=%d" % out_row["end"],
                         "SVLEN=%d" % out_row["svlen"],
                         # CIPOS=-56,20;CIEND=-10,62
                        ])

        yield (out_row["chromosome"], out_row["start"], '.', 'N',
               "<%s>" % out_row["svtype"], '.', '.',
               info, out_row["format"], genotype)


# _____________________________________________________________________________
# THetA

def export_theta(tumor, reference):
    """Convert tumor segments and normal .cnr or reference .cnn to THetA input.

    Follows the THetA segmentation import script but avoid repeating the
    pileups, since we already have the mean depth of coverage in each target
    bin.

    The options for average depth of coverage and read length do not matter
    crucially for proper operation of THetA; increased read counts per bin
    simply increase the confidence of THetA's results.

    THetA2 input format is tabular, with columns:
        ID, chrm, start, end, tumorCount, normalCount

    where chromosome IDs ("chrm") are integers 1 through 24.
    """
    tumor_segs = CNA.read(tumor)
    ref_vals = CNA.read(reference)

    outheader = ["#ID", "chrm", "start", "end", "tumorCount", "normalCount"]
    outrows = []
    # Convert chromosome names to 1-based integer indices
    prev_chrom = None
    chrom_id = 0
    for seg, ref_rows in ref_vals.by_ranges(tumor_segs):
        if seg["chromosome"] != prev_chrom:
            chrom_id += 1
            prev_chrom = seg["chromosome"]
        fields = calculate_theta_fields(seg, ref_rows, chrom_id)
        outrows.append(fields)

    return outheader, outrows


def calculate_theta_fields(seg, ref_rows, chrom_id):
    """Convert a segment's info to a row of THetA input.

    For the normal/reference bin count, take the mean of the bin values within
    each segment so that segments match between tumor and normal.
    """
    # These two scaling factors don't meaningfully affect THetA's calculation
    # unless they're too small
    expect_depth = 100  # Average exome-wide depth of coverage
    read_length = 100
    # Similar number of reads in on-, off-target bins; treat them equally
    segment_size = 1000 * seg["probes"]

    def logratio2count(log2_ratio):
        """Calculate a segment's read count from log2-ratio.

        Math:
            nbases = read_length * read_count
        and
            nbases = segment_size * read_depth
        where
            read_depth = read_depth_ratio * expect_depth

        So:
            read_length * read_count = segment_size * read_depth
            read_count = segment_size * read_depth / read_length
        """
        read_depth = (2 ** log2_ratio) * expect_depth
        read_count = segment_size * read_depth / read_length
        return int(round(read_count))

    tumor_count = logratio2count(seg["log2"])
    ref_count = logratio2count(ref_rows["log2"].mean())
    # e.g. "start_1_93709:end_1_19208166"
    row_id = ("start_%d_%d:end_%d_%d"
              % (chrom_id, seg["start"], chrom_id, seg["end"]))
    return (row_id,       # ID
            chrom_id,     # chrm
            seg["start"], # start
            seg["end"],   # end
            tumor_count,  # tumorCount
            ref_count     # normalCount
           )


# _____________________________________________________________________________

EXPORT_FORMATS = {
    'cdt': fmt_cdt,
    # 'gct': fmt_gct,
    'jtv': fmt_jtv,
    'nexus-basic': export_nexus_basic,
    'nexus-ogt': export_nexus_ogt,
    'seg': export_seg,
    'theta': export_theta,
    'vcf': export_vcf,
}
