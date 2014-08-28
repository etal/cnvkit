#!/usr/bin/env python

"""Extract canonical gene coordinates from UCSC's refFlat.txt.

Usage:
    refFlat2bed.py /path/to/refFlat.txt > genes.bed
    refFlat2bed.py -e /path/to/refFlat.txt > exons.bed
"""
# XXX what about overlapping genes w/ different names?
from __future__ import print_function

from collections import defaultdict

from cnvlib.core import sorter_chrom


def parse_refflat_line(line):
    """Parse one line of refFlat.txt; return relevant info.

    Pair up the exon start and end positions. Add strand direction to each
    chromosome as a "+"/"-" suffix (it will be stripped off later).
    """
    fields = line.split('\t')
    name, refseq, chrom, strand = fields[:4]
    start, end = fields[4:6]
    # start, end = fields[6:8]
    exon_count, exon_starts, exon_ends = fields[8:11]
    exons = zip(map(int, exon_starts.rstrip(',').split(',')),
                map(int, exon_ends.rstrip('\n,').split(',')))
    assert len(exons) == int(exon_count), (
        "Exon count mismatch at %s: file says %s, but counted %d intervals"
        % (name, exon_count, len(exons)))
    return name, refseq, chrom + strand, int(start), int(end), exons


def load_genes(fname):
    genedict = defaultdict(list)
    with open(fname) as genefile:
        for line in genefile:
            name, _refseq, chrom, start, end, _exons = parse_refflat_line(line)
            genedict[name].append((chrom, start, end))
    return genedict


def load_exons(fname):
    exondict = defaultdict(list)
    with open(fname) as genefile:
        for line in genefile:
            name, _refseq, chrom, _start, _end, exons = parse_refflat_line(line)
            for start, end in exons:
                exondict[name].append((chrom, start, end))
    return exondict


def dedupe_regions(regions):
    """Merge genes/exons with the same name and overlapping regions."""
    for name, locs in regions.iteritems():
        locs = list(set(locs))
        if len(locs) == 1:
            chrom, start, end = locs[0]
            # Strip the strand indicator off the end of the chromosome name
            yield (chrom[:-1], start, end, name)
        else:
            # Check locs for overlap; merge if so
            locs.sort()
            chroms, starts, ends = zip(*locs)
            curr_chrom, curr_start, curr_end = chroms[0], starts[0], ends[0]
            for chrom, start, end in zip(chroms[1:], starts[1:], ends[1:]):
                if (chrom != curr_chrom or
                    start > curr_end + 1 or
                    end < curr_start - 1):
                    yield (curr_chrom[:-1], curr_start, curr_end, name)
                    curr_chrom, curr_start, curr_end = chrom, start, end
                else:
                    curr_start = min(start, curr_start)
                    curr_end = max(end, curr_end)
            # Emit the final gene/exon in this group
            yield (curr_chrom[:-1], curr_start, curr_end, name)

def merge_rows(rows):
    """Where rows overlap, merge the regions and gene names."""
    rows = iter(rows)
    prev_chrom, prev_start, prev_end, prev_name = next(rows)

    for chrom, start, end, name in rows:
        # No overlap.
        if (chrom != prev_chrom) or (start >= prev_end):
            yield (prev_chrom, prev_start, prev_end, prev_name)
            # out_row = (chrom, start, end, name)
            prev_chrom, prev_start, prev_end, prev_name = \
                    (chrom, start, end, name)
            continue

        # Some overlap. Adjust prev_ values accordingly.
        # Known: chrom == prev_chrom; start <= prev_end
        assert prev_start <= start, (
            "Botched overlap: %s %s:%s-%s vs. prev. %s %s:%s-%d"
            % (name, chrom, start, end,
               prev_name, prev_chrom, prev_start, prev_end))
        prev_end = max(prev_end, end)
        if name not in prev_name.split('|'):
            prev_name += '|' + name

    # Remainder
    yield (prev_chrom, prev_start, prev_end, prev_name)


def key_genomic_position(row):
    """Turn genomic position into a sort key: (chrom_key, end_posn)

    Input rows are BED-like: (chrom, start, end, name)
    """
    return sorter_chrom(row[0]), row[1]


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('genes', help="refFlat.txt")
    AP.add_argument('-e', '--exons', action='store_true',
                    help="Emit each exon, not just the genes.")
    AP.add_argument('-m', '--merge', action='store_true',
                    help="Merge overlapping regions with different names.")
    args = AP.parse_args()

    if args.exons:
        regions = load_exons(args.genes)
    else:
        regions = load_genes(args.genes)
    out_rows = sorted(dedupe_regions(regions), key=key_genomic_position)
    if args.merge:
        out_rows = merge_rows(out_rows)
    for row in out_rows:
        print('\t'.join(map(str, row)))

