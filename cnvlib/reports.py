"""Supporting functions for the text/tabular-reporting commands.

Namely: breaks, genemetrics.
"""
from __future__ import absolute_import, division
from builtins import str
import collections
import math
import sys

import numpy as np
import pandas as pd

from . import params
from .segmetrics import segment_mean

iteritems = (dict.iteritems if sys.version_info[0] < 3 else dict.items)


# _____________________________________________________________________________
# breaks

def do_breaks(probes, segments, min_probes=1):
    """List the targeted genes in which a copy number breakpoint occurs."""
    intervals = get_gene_intervals(probes)
    bpoints = get_breakpoints(intervals, segments, min_probes)
    return pd.DataFrame.from_records(bpoints,
                                     columns=['gene', 'chromosome',
                                              'location', 'change',
                                              'probes_left', 'probes_right'])


def get_gene_intervals(all_probes, ignore=params.IGNORE_GENE_NAMES):
    """Tally genomic locations of each targeted gene.

    Return a dict of chromosomes to a list of tuples: (gene name, starts, end),
    where gene name is a string, starts is a sorted list of probe start
    positions, and end is the last probe's end position as an integer. (The
    endpoints are redundant since probes are adjacent.)
    """
    ignore += params.ANTITARGET_ALIASES
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
    """Identify segment breaks within the targeted intervals."""
    # TODO use segments.by_ranges(intervals)
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
# genemetrics

def do_genemetrics(cnarr, segments=None, threshold=0.2, min_probes=3,
                skip_low=False, male_reference=False, is_sample_female=None):
    """Identify targeted genes with copy number gain or loss."""
    if is_sample_female is None:
        is_sample_female = cnarr.guess_xx(male_reference=male_reference)
    cnarr = cnarr.shift_xx(male_reference, is_sample_female)
    if segments:
        segments = segments.shift_xx(male_reference, is_sample_female)
        rows = gene_metrics_by_segment(cnarr, segments, threshold, skip_low)
    else:
        rows = gene_metrics_by_gene(cnarr, threshold, skip_low)
    rows = list(rows)
    columns = (rows[0].index if len(rows) else cnarr._required_columns)
    columns = ["gene"] + [col for col in columns if col != "gene"]
    table = pd.DataFrame.from_records(rows).reindex(columns=columns)
    if min_probes and len(table):
        n_probes = (table.segment_probes
                    if 'segment_probes' in table.columns
                    else table.n_bins)
        table = table[n_probes >= min_probes]
    return table


def gene_metrics_by_gene(cnarr, threshold, skip_low=False):
    """Identify genes where average bin copy ratio value exceeds `threshold`.

    NB: Adjust the sample's sex-chromosome log2 values beforehand with shift_xx,
    otherwise all chrX/chrY genes may be reported gained/lost.
    """
    for row in group_by_genes(cnarr, skip_low):
        if abs(row.log2) >= threshold and row.gene:
            yield row


def gene_metrics_by_segment(cnarr, segments, threshold, skip_low=False):
    """Identify genes where segmented copy ratio exceeds `threshold`.

    In the output table, show each segment's weight and probes as segment_weight
    and segment_probes, alongside the gene-level weight and probes.

    NB: Adjust the sample's sex-chromosome log2 values beforehand with shift_xx,
    otherwise all chrX/chrY genes may be reported gained/lost.
    """
    extra_cols = [col for col in segments.data.columns
                  if col not in cnarr.data.columns
                  and col not in ('depth', 'probes', 'weight')]
    for colname in extra_cols:
        cnarr[colname] = np.nan
    for segment, subprobes in cnarr.by_ranges(segments):
        if abs(segment.log2) >= threshold:
            for row in group_by_genes(subprobes, skip_low):
                row["log2"] = segment.log2
                if hasattr(segment, 'weight'):
                    row['segment_weight'] = segment.weight
                if hasattr(segment, 'probes'):
                    row['segment_probes'] = segment.probes
                for colname in extra_cols:
                    row[colname] = getattr(segment, colname)
                yield row


# ENH consolidate with CNA.squash_genes
def group_by_genes(cnarr, skip_low):
    """Group probe and coverage data by gene.

    Return an iterable of genes, in chromosomal order, associated with their
    location and coverages:

        [(gene, chrom, start, end, [coverages]), ...]
    """
    ignore = ('', np.nan) + params.ANTITARGET_ALIASES
    for gene, rows in cnarr.by_gene():
        if not rows or gene in ignore:
            continue
        segmean = segment_mean(rows, skip_low)
        if segmean is None:
            continue
        outrow = rows[0].copy()
        outrow["end"] = rows.end.iat[-1]
        outrow["gene"] = gene
        outrow["log2"] = segmean
        outrow["n_bins"] = len(rows)
        if "weight" in rows:
            outrow["weight"] = rows["weight"].sum()
            if "depth" in rows:
                outrow["depth"] = np.average(rows["depth"],
                                             weights=rows["weight"])
        elif "depth" in rows:
            outrow["depth"] = rows["depth"].mean()
        yield outrow
