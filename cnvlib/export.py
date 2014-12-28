"""Export CNVkit objects and files to other formats."""
from __future__ import absolute_import, division, print_function
import collections
import sys

from Bio._py3k import map, range, zip

from . import core
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
# freebayes

def export_freebayes(sample_fname, args):
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

    segs = CNA.read(sample_fname)
    print(args.gender)
    is_sample_female = core.guess_xx(segs, args.male_normal, verbose=False)
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

    bedrows = segments2bed(segs, args.name or segs.sample_id,
                           args.ploidy, args.purity, args.male_normal,
                           is_sample_female)
    return None, list(bedrows)


def segments2bed(segments, sample_name, ploidy, purity, is_reference_male,
                 is_sample_female):
    for row in segments:
        ref_copies, expect_copies = reference_expect_copies(
            row["chromosome"], ploidy, is_sample_female, is_reference_male)
        ncopies = log2_ratio_to_integer(
            row["coverage"], ref_copies, expect_copies, purity)
        # Ignore regions of neutral copy number
        if ncopies != ploidy:
            yield (row["chromosome"], # reference sequence
                   row["start"], # start (0-indexed)
                   row["end"], # end
                   sample_name, # sample name
                   ncopies) # copy number


def reference_expect_copies(chrom, ploidy, is_sample_female, is_reference_male):
    """Determine the number copies of a chromosome expected and in reference.

    For sex chromosomes, these values may not be the same ploidy as the
    autosomes.

    Return a pair: number of copies in the reference and expected in the sample.
    """
    chrom = chrom.lower()

    if chrom in ["chrx", "x"]:
        ref_copies = (ploidy // 2 if is_reference_male else ploidy)
        exp_copies = (ploidy if is_sample_female else ploidy // 2)
    elif chrom in ["chry", "y"]:
        ref_copies = ploidy // 2
        exp_copies = (0 if is_sample_female else ploidy // 2)
    else:
        ref_copies = exp_copies = ploidy

    return ref_copies, exp_copies


def log2_ratio_to_integer(log2_ratio, ref_copies, expect_copies, purity=None):
    """Transform a log2 ratio value to integer.

    Math:

        log2_ratio = log2(ncopies / ploidy)
        2^log2_ratio = ncopies / ploidy
        ncopies = ploidy * 2^log2_ratio

    With rescaling for purity:

        let v = log2 ratio value, p = tumor purity,
            r = reference ploidy, x = expected ploidy;
        v = log_2(p*n/r + (1-p)*x/r)
        2^v = p*n/r + (1-p)*x/r
        n*p/r = 2^v - (1-p)*x/r
        n = (r*2^v - x*(1-p)) / p

    If purity adjustment is skipped (p=1; e.g. if THetA was run beforehand):

        n = r*2^v
    """
    EPSILON = 1e-7

    if purity and purity < 1.0:
        ncopies = (ref_copies * 2**log2_ratio - expect_copies * (1 - purity)
                  ) / purity
    else:
        ncopies = ref_copies * 2 ** log2_ratio
        # Convention: encode log2(0 copies) as a half-copy (-2 for diploid)
        if ncopies <= .5 + EPSILON:
            return 0
    return max(0, int(round(ncopies)))


# _____________________________________________________________________________
# theta

def export_theta(sample_fname):
    pass


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
