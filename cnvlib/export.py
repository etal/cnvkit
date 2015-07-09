"""Export CNVkit objects and files to other formats."""
from __future__ import absolute_import, division, print_function

import collections
import math
import sys

from Bio._py3k import map, range, zip

from . import call, core
from .cnarray import row2label, CopyNumArray as CNA

ProbeInfo = collections.namedtuple('ProbeInfo', 'label chrom start end gene')

def merge_samples(filenames):
    """Merge probe values from multiple samples into a 2D table (of sorts).

    Input:
        dict of {sample ID: (probes, values)}
    Output:
        list-of-tuples: (probe, log2 coverages...)
    """
    handles = []
    datastreams = []  # e.g. [list-of-pairs, list-of-pairs, ...]
    for fname in filenames:
        handle = open(fname)
        handles.append(handle)
        data = core.parse_tsv(handle)
        datastreams.append(data)
    # Emit the individual rows merged across samples, one per probe
    for rows in zip(*datastreams):
        yield merge_rows(rows)
    # Clean up
    for handle in handles:
        handle.close()


def merge_rows(rows):
    """Combine equivalent rows of coverage data across multiple samples.

    Check that probe info matches across all samples, then merge the log2
    coverage values.

    Input: a list of individual rows corresponding to the same probes from
    different coverage files.
    Output: a list starting with the single common Probe object, followed by the
    log2 coverage values from each sample, in order.
    """
    probe_infos, coverages = zip(*map(row_to_probe_coverage, rows))
    probe_info = core.check_unique(probe_infos, "probe Name")
    combined_row = [probe_info] + list(coverages)
    return combined_row


def row_to_probe_coverage(row):
    """Repack a parsed row into a ProbeInfo instance and coverage value."""
    chrom, start, end, gene, coverage = row[:5]
    label = "%s:%s-%s:%s" % (chrom, start, end, gene)
    probe_info = ProbeInfo(label, chrom, int(start), int(end), gene)
    return probe_info, float(coverage)


# Supported formats:

def fmt_cdt(sample_ids, rows):
    """Format as CDT."""
    outheader = ['GID', 'CLID', 'NAME', 'GWEIGHT'] + sample_ids
    header2 = ['AID', '', '', '']
    header2.extend(['ARRY' + str(i).zfill(3) + 'X'
                    for i in range(len(sample_ids))])
    outrows = [header2]
    for i, row in enumerate(rows):
        probe, values = row[0], row[1:]
        outrow = ['GENE%dX' % i, 'IMAGE:%d' % i, probe.label, 1] # or probe.gene?
        outrow.extend(values)
        outrows.append(outrow)
    return outheader, outrows


# TODO
def fmt_gct(sample_ids, rows):
    return NotImplemented


def fmt_jtv(sample_ids, rows):
    """Format for Java TreeView."""
    outheader = ["CloneID", "Name"] + sample_ids
    outrows = [["IMAGE:", row[0].label] + row[1:] for row in rows]
    return outheader, outrows


# TODO
def fmt_multi(sample_ids, rows):
    return NotImplemented


# TODO
def fmt_vcf(sample_ids, rows):
    return NotImplemented


# Special cases

def export_nexus_basic(sample_fname):
    """Biodiscovery Nexus Copy Number "basic" format.

    Only represents one sample per file.
    """
    cnarr = CNA.read(sample_fname)
    outheader = ['probe', 'chromosome', 'start', 'end', 'gene', 'log2']
    outrows = [[row2label(row)] + list(row)[:5] for row in cnarr.data]
    return outheader, outrows


def export_seg(sample_fnames):
    """SEG format for copy number segments.

    Segment breakpoints are not the same across samples, so samples are listed
    in serial with the sample ID as the left column.
    """
    outrows = []
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

        if 'probes' in segments.data.dtype.fields:
            outheader = ["ID", "Chromosome", "Start", "End", "NumProbes", "Mean"]
            def row2out(row):
                return (segments.sample_id, chrom_ids[row['chromosome']],
                        row['start'], row['end'], row['probes'],
                        row['coverage'])
        else:
            outheader = ["ID", "Chromosome", "Start", "End", "Mean"]
            def row2out(row):
                return (segments.sample_id, chrom_ids[row['chromosome']],
                        row['start'], row['end'], row['coverage'])
        outrows.extend(row2out(row) for row in segments)
    return outheader, outrows


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

def export_bed(sample_fnames, args):
    """Export to BED format.

    For each region in each sample which does not have neutral copy number
    (equal to 2 or the value set by --ploidy), the columns are:

        - reference sequence name
        - start (0-indexed)
        - end
        - sample name or given label
        - integer copy number
    """
    bed_rows = []
    for fname in sample_fnames:
        segs = CNA.read(fname)
        rows = segments2bed(segs, args.sample_id or segs.sample_id, args.ploidy,
                            args.male_reference)
        bed_rows.extend(rows)
    return None, bed_rows


def segments2bed(segments, label, ploidy, is_reference_male):
    """Convert a copy number array to a BED-like format."""
    absolutes = call.absolute_pure(segments, ploidy, is_reference_male)
    for row, abs_val in zip(segments, absolutes):
        ncopies = int(round(abs_val))
        # Ignore regions of neutral copy number
        if ncopies != ploidy:
            yield (row["chromosome"], row["start"], row["end"], label, ncopies)


# _____________________________________________________________________________
# freebayes

def export_freebayes(sample_fnames, args):
    """Export to FreeBayes --cnv-map format.

    Which is BED-like, for each region in each sample which does not have
    neutral copy number (equal to 2 or the value set by --ploidy), with columns:

        - reference sequence
        - start (0-indexed)
        - end
        - sample name
        - copy number
    """
    if args.purity and not 0.0 < args.purity <= 1.0:
        raise RuntimeError("Purity must be between 0 and 1.")

    bed_rows = []
    for fname in sample_fnames:
        segs = CNA.read(fname)
        is_sample_female = core.guess_xx(segs, args.male_reference,
                                         verbose=False)
        if args.gender:
            is_sample_female_given = (args.gender in ["f", "female"])
            if is_sample_female != is_sample_female_given:
                print("Sample gender specified as", args.gender,
                      "but chrX copy number looks like",
                      "female" if is_sample_female else "male",
                      file=sys.stderr)
                is_sample_female = is_sample_female_given
        print("Treating sample gender as",
              "female" if is_sample_female else "male",
              file=sys.stderr)
        rows = segments2freebayes(segs, args.sample_id or segs.sample_id,
                                  args.ploidy, args.purity, args.male_reference,
                                  is_sample_female)
        bed_rows.extend(rows)
    return None, bed_rows


def segments2freebayes(segments, sample_name, ploidy, purity, is_reference_male,
                       is_sample_female):
    """Convert a copy number array to a BED-like format."""
    absolutes = call.absolute_clonal(segments, ploidy, purity,
                                     is_reference_male, is_sample_female)
    for row, abs_val in zip(segments, absolutes):
        ncopies = int(round(abs_val))
        # Ignore regions of neutral copy number
        if ncopies != ploidy:
            yield (row["chromosome"], # reference sequence
                   row["start"], # start (0-indexed)
                   row["end"], # end
                   sample_name, # sample name
                   ncopies) # copy number


# _____________________________________________________________________________
# theta

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
    for seg, ref_rows in ref_vals.by_segment(tumor_segs):
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

    tumor_count = logratio2count(seg["coverage"])
    ref_count = logratio2count(ref_rows["coverage"].mean())
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
    # 'vcf': fmt_vcf,
    'nexus-basic': export_nexus_basic,
    # 'nexus-multi1': fmt_multi,
    'seg': export_seg,
    'freebayes': export_freebayes,
    'theta': export_theta,
}
