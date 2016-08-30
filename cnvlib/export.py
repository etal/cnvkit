"""Export CNVkit objects and files to other formats."""
from __future__ import absolute_import, division, print_function
from builtins import map, range, str, zip

import collections
import logging
import time

import numpy as np
import pandas as pd

from . import call, core, params, tabio
from ._version import __version__


def merge_samples(filenames):
    """Merge probe values from multiple samples into a 2D table (of sorts).

    Input:
        dict of {sample ID: (probes, values)}
    Output:
        list-of-tuples: (probe, log2 coverages...)
    """
    def label_with_gene(cnarr):
        row2label = lambda row: "{}:{}-{}:{}".format(
            row.chromosome, row.start, row.end, row.gene)
        return cnarr.data.apply(row2label, axis=1)

    if not filenames:
        return []
    first_cnarr = tabio.read_cna(filenames[0])
    out_table = first_cnarr.data.loc[:, ["chromosome", "start", "end", "gene"]]
    out_table["label"] = label_with_gene(first_cnarr)
    out_table[first_cnarr.sample_id] = first_cnarr["log2"]
    for fname in filenames[1:]:
        cnarr = tabio.read_cna(fname)
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
    cnarr = tabio.read_cna(sample_fname)
    out_table = cnarr.data.loc[:, ['chromosome', 'start', 'end', 'gene', 'log2']]
    out_table['probe'] = cnarr.labels()
    return out_table


def export_nexus_ogt(sample_fname, vcf_fname, sample_id,
                     min_depth=20, min_weight=0.0):
    """Biodiscovery Nexus Copy Number "Custom-OGT" format.

    To create the b-allele frequencies column, alterate allele frequencies from
    the VCF are aligned to the .cnr file bins.  Bins that contain no variants
    are left blank; if a bin contains multiple variants, then the frequencies
    are all "mirrored" to be above or below .5 (majority rules), then the median
    of those values is taken.
    """
    cnarr = tabio.read_cna(sample_fname)
    if min_weight and "weight" in cnarr:
        mask_low_weight = (cnarr["weight"] < min_weight)
        logging.info("Dropping %d bins with weight below %f",
                     mask_low_weight.sum(), min_weight)
        cnarr.data = cnarr.data[~mask_low_weight]
    varr = tabio.read(vcf_fname, "vcf",
                      sample_id=sample_id or cnarr.sample_id,
                      min_depth=min_depth, skip_somatic=True)
    bafs = varr.heterozygous().baf_by_ranges(cnarr)
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
        segments = tabio.read_cna(fname)
        if chrom_ids is None:
            # Create & store
            chrom_ids = create_chrom_ids(segments)
        else:
            # Verify
            core.assert_equal("Segment chromosome names differ",
                              previous=list(chrom_ids.keys()),
                              current=list(create_chrom_ids(segments).keys()))
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

def export_bed(segments, ploidy, is_reference_male, is_sample_female,
               label, show):
    """Convert a copy number array to a BED-like DataFrame.

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
    out = segments.data.loc[:, ["chromosome", "start", "end"]]
    out["label"] = label
    out["ncopies"] = (segments["cn"] if "cn" in segments
                      else np.rint(call.absolute_pure(segments, ploidy,
                                                      is_reference_male)))
    if show == "ploidy":
        # Skip regions of default ploidy
        out = out[out["ncopies"] != ploidy]
    elif show == "variant":
        # Skip regions of non-neutral copy number
        exp_copies = call.absolute_expect(segments, ploidy, is_sample_female)
        out = out[out["ncopies"] != exp_copies]
    return out


# _____________________________________________________________________________
# VCF

VCF_HEADER = """\
##fileformat=VCFv4.0
##fileDate={date}
##source=CNVkit v{version}
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=FOLD_CHANGE,Number=1,Type=Float,Description="Fold change">
##INFO=<ID=FOLD_CHANGE_LOG,Number=1,Type=Float,Description="Log fold change">
##INFO=<ID=PROBES,Number=1,Type=Integer,Description="Number of probes in CNV">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
""".format(date=time.strftime("%Y%m%d"), version=__version__)
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
    if "cn" in segments:
        out_dframe["ncopies"] = segments["cn"]
        abs_expect = call.absolute_expect(segments, ploidy, is_sample_female)

    else:
        abs_dframe = call.absolute_dataframe(segments, ploidy, 1.0,
                                             is_reference_male,
                                             is_sample_female)
        out_dframe["ncopies"] = np.rint(abs_dframe["absolute"])
        abs_expect = abs_dframe["expect"]
    idx_losses = (out_dframe["ncopies"] < abs_expect)

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
    for out_row, abs_exp in zip(out_dframe.itertuples(index=False), abs_expect):
        if (out_row.ncopies == abs_exp or
            # Survive files from buggy v0.7.1 (#53)
            not str(out_row.probes).isdigit()):
            # Skip regions of neutral copy number
            continue  # or "CNV" for subclonal?

        if out_row.ncopies > abs_exp:
            genotype = "0/1:0:%d:%d" % (out_row.ncopies, out_row.probes)
        elif out_row.ncopies < abs_exp:
            # TODO XXX handle non-diploid ploidies, haploid chroms
            if out_row.ncopies == 0:
                # Complete deletion, 0 copies
                gt = "1/1"
            else:
                # Single copy deletion
                gt = "0/1"
            genotype = "%s:%d" % (gt, out_row.probes)

        info = ";".join(["IMPRECISE",
                         "SVTYPE=%s" % out_row.svtype,
                         "END=%d" % out_row.end,
                         "SVLEN=%d" % out_row.svlen,
                         "FOLD_CHANGE=%f" % 2.0 ** out_row.log2,
                         "FOLD_CHANGE_LOG=%f" % out_row.log2,
                         "PROBES=%d" % out_row.probes
                         # CIPOS=-56,20;CIEND=-10,62
                        ])

        yield (out_row.chromosome, out_row.start, '.', 'N',
               "<%s>" % out_row.svtype, '.', '.',
               info, out_row.format, genotype)


# _____________________________________________________________________________
# THetA

def export_theta(tumor_segs, normal_cn):
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
    out_columns = ["#ID", "chrm", "start", "end", "tumorCount", "normalCount"]
    if not tumor_segs:
        return pd.DataFrame(columns=out_columns)

    # Drop any chromosomes that are not integer or XY
    # THetA hard-codes X & Y, so we can, too
    xy_names = (["chrX", "chrY"]
                if tumor_segs.chromosome.iat[0].startswith('chr')
                else ["X", "Y"])
    tumor_segs = tumor_segs.autosomes(also=xy_names)
    if normal_cn:
        normal_cn = normal_cn.autosomes(also=xy_names)

    table = tumor_segs.data.loc[:, ["start", "end"]]

    # Convert chromosome names to 1-based integer indices
    chr2idx = {c: i+1
            for i, c in enumerate(tumor_segs.chromosome.drop_duplicates())}
    table["chrm"] = tumor_segs.chromosome.replace(chr2idx)
    # Unique string identifier for each row, e.g. "start_1_93709:end_1_19208166"
    table["#ID"] = ["start_%d_%d:end_%d_%d"
                    % (row.chrm, row.start, row.chrm, row.end)
                    for row in table.itertuples(index=False)]

    # Calculate/estimate per-segment read counts in tumor and normal samples
    ref_means, nbins = ref_means_nbins(tumor_segs, normal_cn)
    table["tumorCount"] = theta_read_counts(tumor_segs.log2, nbins)
    table["normalCount"] = theta_read_counts(ref_means, nbins)
    return table[out_columns]


def ref_means_nbins(tumor_segs, normal_cn):
    """Extract segments' reference mean log2 values and probe counts.

    Code paths::

        wt_mdn  wt_old  probes  norm -> norm, nbins
        +       *       *       -       0,  wt_mdn
        -       +       +       -       0,  wt_old * probes
        -       +       -       -       0,  wt_old * size?
        -       -       +       -       0,  probes
        -       -       -       -       0,  size?

        +       -       +       +       norm, probes
        +       -       -       +       norm, bin counts
        -       +       +       +       norm, probes
        -       +       -       +       norm, bin counts
        -       -       +       +       norm, probes
        -       -       -       +       norm, bin counts
    """
    if normal_cn:
        subarrs = [subarr for _seg, subarr in normal_cn.by_ranges(tumor_segs)]
        # For the normal/reference bin count, take the mean of the bin values
        # within each segment so that segments match between tumor and normal.
        # ENH: weighted mean, like gainloss
        ref_means = np.asarray([s.log2.mean() for s in subarrs])
        if "probes" in tumor_segs:
            nbins = tumor_segs["probes"]
        else:
            nbins = np.asarray([len(s) for s in subarrs])
    else:
        # Assume neutral reference log2 across all segments
        # (weights already account for reference log2)
        ref_means = np.zeros(len(tumor_segs))
        if "weight" in tumor_segs and (tumor_segs["weight"] > 1.0).any():
            # Segment weights are already multiplied by probe counts
            nbins = tumor_segs["weight"]
            # Rescale to average 1.0 (we'll do the same below, too)
            nbins /= nbins.max() / nbins.mean()
        else:
            if "probes" in tumor_segs:
                nbins = tumor_segs["probes"]
            else:
                logging.warn("No probe counts in tumor segments file and no "
                             "normal reference given; guessing normal "
                             "read-counts-per-segment from segment sizes")
                sizes = tumor_segs.end - tumor_segs.start
                nbins = sizes / sizes.mean()
            if "weight" in tumor_segs:
                # Already checked -- these are old-style weights that were not
                # multiplied by `probes` originally
                nbins *= tumor_segs["weight"] / tumor_segs["weight"].mean()
    return ref_means, nbins


def theta_read_counts(log2_ratio, nbins,
                      # These two scaling factors don't meaningfully affect
                      # THetA's calculation unless they're too small
                      avg_depth=500, avg_bin_width=200):
    """Calculate segments' read counts from log2-ratios.

    Math:
        nbases = read_length * read_count
    and
        nbases = bin_width * read_depth
    where
        read_depth = read_depth_ratio * avg_depth

    So:
        read_length * read_count = bin_width * read_depth
        read_count = bin_width * read_depth / read_length
    """
    read_depth = (2 ** log2_ratio) * avg_depth
    read_count = nbins * avg_bin_width * read_depth / params.READ_LEN
    return read_count.round().astype('int')


def export_theta_snps(varr):
    """Generate THetA's SNP per-allele read count "formatted.txt" files."""
    # Drop any chromosomes that are not integer or XY
    varr = varr.autosomes(also=(["chrX", "chrY"]
                       if varr.chromosome.iat[0].startswith("chr")
                       else ["X", "Y"]))
    # Skip indels
    varr = varr[(varr["ref"].str.len() == 1) & (varr["alt"].str.len() == 1)]
    # Drop rows with any NaN
    varr.data.dropna(subset=["depth", "alt_count"], inplace=True)
    if "n_depth" in varr and "n_alt_count" in varr:
        varr.data.dropna(subset=["n_depth", "alt_count"], inplace=True)
    # Avoid weird situation I've seen on alt contigs
    varr = varr[varr["depth"] >= varr["alt_count"]]
    # Reformat for THetA2
    for depth_key, alt_key in (("depth", "alt_count"),
                               ("n_depth", "n_alt_count")):
        # Extract the SNP allele counts that THetA2 uses
        table = varr.data.loc[:, ("chromosome", "start", depth_key, alt_key)]
        table["ref_depth"] = (table[depth_key] - table[alt_key]).astype("int")
        table[alt_key] = table[alt_key].astype("int")
        table = table.loc[:, ("chromosome", "start", "ref_depth", alt_key)]
        table.columns = ["#Chrm", "Pos", "Ref_Allele", "Mut_Allele"]
        yield table


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
