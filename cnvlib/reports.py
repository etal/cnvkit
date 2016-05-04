"""Supporting functions for the text/tabular-reporting commands.

Namely: breaks, gainloss.
"""
from __future__ import absolute_import, division
from builtins import str
import collections
import math
import sys

from . import metrics, params

iteritems = (dict.iteritems if sys.version_info[0] < 3 else dict.items)


# _____________________________________________________________________________
# breaks

def get_gene_intervals(all_probes, ignore=params.IGNORE_GENE_NAMES):
    """Tally genomic locations of each targeted gene.

    Return a dict of chromosomes to a list of tuples: (gene name, start, end).
    """
    ignore += ("Background",)
    # Tally the start & end points for each targeted gene; group by chromosome
    gene_probes = collections.defaultdict(lambda: collections.defaultdict(list))
    for row in all_probes:
        gname = str(row.gene)
        if gname not in ignore:
            gene_probes[row.chromosome][gname].append(row)
    # Condense into a single interval for each gene
    intervals = collections.defaultdict(list)
    for chrom, gp in iteritems(gene_probes):
        for gene, probes in iteritems(gp):
            starts = sorted(row.start for row in probes)
            end = max(row.end for row in probes)
            intervals[chrom].append((gene, starts, end))
        intervals[chrom].sort(key=lambda gse: gse[1])
    return intervals


def get_breakpoints(intervals, segments, min_probes):
    """Identify CBS segment breaks within the targeted intervals."""
    breakpoints = []
    for i, curr_row in enumerate(segments[:-1]):
        curr_chrom = curr_row.chromosome
        curr_end = curr_row.end
        next_row = segments[i + 1]
        # Skip if this segment is the last (or only) one on this chromosome
        if next_row.chromosome != curr_chrom:
            continue
        for gname, gstarts, gend in intervals[curr_chrom]:
            if gstarts[0] < curr_end < gend:
                probes_left = sum(s < curr_end for s in gstarts)
                probes_right = sum(s >= curr_end for s in gstarts)
                if probes_left >= min_probes and probes_right >= min_probes:
                    breakpoints.append(
                        (gname, curr_chrom, int(math.ceil(curr_end)),
                         next_row.log2 - curr_row.log2,
                         probes_left, probes_right))
    breakpoints.sort(key=lambda row: (min(row[4], row[5]), abs(row[3])),
                     reverse=True)
    return breakpoints


# _____________________________________________________________________________
# gainloss

def gainloss_by_gene(probes, threshold, skip_low=False):
    """Identify genes where average bin copy ratio value exceeds `threshold`.

    NB: Must shift sex-chromosome values beforehand with shift_xx,
    otherwise all chrX/chrY genes may be reported gained/lost.
    """
    for gene, chrom, start, end, coverage, nprobes in group_by_genes(probes,
                                                                     skip_low):
        if abs(coverage) >= threshold:
            yield (gene, chrom, start, end, coverage, nprobes)


def gainloss_by_segment(probes, segments, threshold, skip_low=False):
    """Identify genes where segmented copy ratio exceeds `threshold`.

    NB: Must shift sex-chromosome values beforehand with shift_xx,
    otherwise all chrX/chrY genes may be reported gained/lost.
    """
    for segment, subprobes in probes.by_ranges(segments):
        if abs(segment.log2) >= threshold:
            for (gene, chrom, start, end, _coverage, nprobes
                ) in group_by_genes(subprobes, skip_low):
                yield (gene, chrom, start, end, segment.log2, nprobes)


# TODO consolidate with CNA.squash_genes
def group_by_genes(probes, skip_low):
    """Group probe and coverage data by gene.

    Return an iterable of genes, in chromosomal order, associated with their
    location and coverages:

        [(gene, chrom, start, end, [coverages]), ...]
    """
    for gene, rows in probes.by_gene():
        if gene == 'Background':
            continue
        chrom = rows[0, 'chromosome']
        start = rows[0, 'start']
        end = rows[len(rows)-1, 'end']
        segmean = metrics.segment_mean(rows)
        if segmean is not None:
            nprobes = len(rows)
            yield gene, chrom, start, end, segmean, nprobes
