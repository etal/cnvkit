"""Export CNVkit objects and files to other formats."""
from __future__ import absolute_import, division
import collections

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


EXPORT_FORMATS = {
    'cdt': fmt_cdt,
    # 'gct': fmt_gct,
    'jtv': fmt_jtv,
    # 'vcf': fmt_vcf,
    'nexus-basic': export_nexus_basic,
    # 'nexus-multi1': fmt_multi,
    'seg': export_seg,
}
