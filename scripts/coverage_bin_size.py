#!/usr/bin/env python
"""Estimate on- and off-target coverages to calculate good average bin sizes."""
from __future__ import absolute_import, division, print_function
import argparse
import logging
from io import BytesIO

import numpy as np
import pandas as pd
import pysam

from cnvlib import samutil, tabio
from cnvlib.descriptives import weighted_median
from cnvlib.genome.gary import GenomicArray as GA

logging.basicConfig(level=logging.WARN, format="%(message)s")

# Desired average number of bases per bin
# E.g.
#  reads-per-bin * read-length == binsize * depth = 50kbp
#  binsize = 50kbp / depth
#  depth = 500 -> binsize = 50,000 / 500  == 100
BP_PER_BIN = 50000.0


def main(args):
    # Validate
    if args.method in ('hybrid', 'amplicon') and not args.targets:
        raise RuntimeError("Sequencing method %r requires targets", args.method)
    elif args.method == 'wgs' and args.targets:
        logging.warn("Targets will be ignored: %s", args.targets)
    if args.method == 'amplicon' and args.access:
        logging.warn("Sequencing-accessible regions will be ignored: %s",
                     args.access)
    samutil.ensure_bam_index(args.bam)
    # Prep
    rc_table = samutil.bam_read_counts(args.bam, drop_unmapped=True)
    #  bam = pysam.Samfile(args.bam)
    read_len = samutil.get_read_length(args.bam)
    logging.warn("Estimated read length %s", read_len)
    access = read_regions(args.access)
    targets = read_regions(args.targets)
    # Dispatch
    if args.method == 'amplicon':
        amplicon(rc_table, read_len, args.bam, targets)
    elif args.method == 'hybrid':
        hybrid(rc_table, read_len, args.bam, targets, access)
    elif args.method == 'wgs':
        wgs(rc_table, read_len, access)


def amplicon(rc_table, read_len, bam_fname, targets):
    """Targeted amplicon sequencing.

    amplicon: -t
        no antitargets or -g; all reads are in targets
        chrom mean cvg = read count / sum(target sizes)
    """
    rc_table = update_chrom_length(rc_table, targets)
    depth = average_depth(rc_table, read_len)
    print("          ", "Depth", "Bin size", sep='\t')
    print("On-target:", depth, BP_PER_BIN / depth, sep='\t')
    # DBG: Compare to sample_region_depths
    #  print("\nALTERNATIVELY --")
    #  depth = sample_region_depths(bam_fname, targets)
    #  print("          ", "Depth", "Bin size", sep='\t')
    #  print("On-target:", depth, BP_PER_BIN / depth, sep='\t')
    #  print("\nALTERNATIVELY --")
    #  depth = sample_region_cov(bam_fname, targets)
    #  print("          ", "Depth", "Bin size", sep='\t')
    #  print("On-target:", depth, BP_PER_BIN / depth, sep='\t')


def hybrid(rc_table, read_len, bam_fname, targets, access=None):
    """Hybrid capture sequencing.

    hybrid: -t, (-g)
        w/o -g, access is (0, chromosome lengths)
        antitargets: access.subtract(targets)
        target size, antitarget size: sum of region sizes
        target/antitarget mean cvg: read count / sizes
            XXX how many of those reads are in targets?
                randomly sample 100 targets & count reads?
                    (in middle 50% of sizes)
    """
    # Only examine chromosomes present in all 2-3 input datasets

    ok_chroms = shared_chroms([rc_table, targets] +
                              [access] if access is not None else [])
    rc_table = rc_table[rc_table.chromosome.isin(ok_chroms)]
    targets = targets[targets.chromosome.isin(ok_chroms)]

    # Identify off-target regions
    if access is None:
        access = idxstats2ga(rc_table)
    else:
        access = access[access.chromosome.isin(ok_chroms)]
    antitargets = access.subtract(targets)
    anti_table = update_chrom_length(rc_table, antitargets)
    # Deal with targets
    target_depth = sample_region_cov(bam_fname, targets)
    # Account for captured reads
    target_length = region_size_by_chrom(targets)['length']
    target_reads = (target_length * target_depth / read_len).values
    anti_table = anti_table.assign(mapped=anti_table.mapped - target_reads)
    anti_depth = average_depth(anti_table, read_len)
    print("           ", "Depth", "Bin size", sep='\t')
    print("On-target: ", target_depth, BP_PER_BIN / target_depth, sep='\t')
    print("Off-target:", anti_depth, BP_PER_BIN / anti_depth, sep='\t')
    # DBG
    #  print("ALTERNATIVELY --")
    #  anti_length = anti_table['length']
    #  anti_depth = sample_region_cov(bam_fname, antitargets)
    #  anti_reads = (anti_length * anti_depth / read_len).values
    #  tgt_table = update_chrom_length(rc_table, targets)
    #  tgt_table = tgt_table.assign(mapped=tgt_table.mapped - anti_reads)
    #  tgt_depth = average_depth(tgt_table, read_len)
    #  print("           ", "Depth", "Bin size", sep='\t')
    #  print("On-target: ", tgt_depth, BP_PER_BIN / tgt_depth, sep='\t')
    #  print("Off-target:", anti_depth, BP_PER_BIN / anti_depth, sep='\t')


def wgs(rc_table, read_len, access=None):
    """Whole genome sequencing.

    wgs: (-g)
        chrom mean cvg = read count / sum(access sizes)
    """
    if access is not None and len(access):
        rc_table = update_chrom_length(rc_table, access)
    depth = average_depth(rc_table, read_len)
    print("       ", "Depth", "Bin size", sep='\t')
    print("Genome:", depth, BP_PER_BIN / depth, sep='\t')


# ---

def read_regions(bed_fname):
    if bed_fname:
        regions = tabio.read_auto(bed_fname)
        if len(regions):
            return regions
        else:
            logging.warn("No regions to estimate depth from %s",
                        regions.meta.get('filename', ''))
    return None


def update_chrom_length(table, regions):
    if regions is not None and len(regions):
        chrom_sizes = region_size_by_chrom(regions)
        table = table.merge(chrom_sizes, on='chromosome', how='inner')
        table['length'] = table['length_y'] # ?
    return table


def region_size_by_chrom(regions):
    chromgroups = regions.data.groupby('chromosome', sort=False)
    # sizes = chromgroups.apply(total_region_size) # XXX
    sizes = [total_region_size(g) for _key, g in chromgroups]
    return pd.DataFrame({'chromosome': regions.chromosome.drop_duplicates(),
                         'length': sizes})


def total_region_size(regions):
    """Aggregate area of all genomic ranges in `regions`."""
    return (regions.end - regions.start).sum()


def average_depth(table, read_length=100):
    """Estimate the average read depth across the genome.

    Returns
    -------
    float
        Median of the per-chromosome mean read depths, weighted by chromosome
        size.
    """
    mean_depths = read_length * table.mapped / table.length
    return weighted_median(mean_depths, table.length)


def sample_region_depths(bam_fname, regions, max_num=100):
    """."""
    bam = pysam.Samfile(bam_fname, 'rb')
    sizes = regions.end - regions.start
    lo_size, hi_size = np.percentile(sizes, [25, 75])
    midsize_regions = regions.data[(sizes >= lo_size) & (sizes <= hi_size)]
    if len(midsize_regions) > max_num:
        midsize_regions = midsize_regions.sample(max_num)
    counter = 0
    for reg in midsize_regions.itertuples(index=False):
        counter += bam.count(reg.chromosome, reg.start, reg.end)
    return counter / len(midsize_regions)


def sample_region_cov(bam_fname, regions, max_num=100):
    """."""
    sizes = regions.end - regions.start
    lo_size, hi_size = np.percentile(sizes, [25, 75])
    midsize_regions = regions.data[(sizes >= lo_size) & (sizes <= hi_size)]
    if len(midsize_regions) > max_num:
        midsize_regions = midsize_regions.sample(max_num)
    # XXX
    import tempfile
    with tempfile.NamedTemporaryFile(suffix='.bed', # mode='w+t'
                                    ) as f:
        #  f.name
        tabio.write(regions.as_dataframe(midsize_regions), f, 'bed4')
        f.flush()
        #  cnvlib.coverage.bedcov(f.name, bam_fname, 0)
        raw = pysam.bedcov(f.name, bam_fname)
        table = pd.read_table(BytesIO(raw),
                              names=['chrom', 'start', 'end', 'gene', 'bp'])
    assert len(table), table
    table = table[(table.end > table.start) & (table.bp > 0)]  # Safety
    assert len(table), table
    depths = table.bp / (table.end - table.start)
    avg_depth = depths.mean()
    assert avg_depth > 0, table
    return avg_depth


def shared_chroms(tables):
    """Intersection of DataFrame .chromosome values."""
    chroms = tables[0].chromosome.drop_duplicates()
    for tab in tables[1:]:
        new_chroms = tab.chromosome.drop_duplicates()
        chroms = chroms[chroms.isin(new_chroms)]
    return chroms


def idxstats2ga(table):
    return GA(table.assign(start=0, end=table.length)
              .loc[:, ('chromosome', 'start', 'end')])



if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('bam', #nargs='+',
                    help="""Sample BAM file to test for target coverage""")
    AP.add_argument('-m', '--method',
                    choices=('hybrid', 'amplicon', 'wgs'), default='hybrid',
                    help="""Sequencing protocol: hybridization capture
                    ('hybrid'), targeted amplicon sequencing ('amplicon'), or
                    whole genome sequencing ('wgs'). Determines whether and how
                    to use antitarget bins. [Default: %(default)s]""")
    AP.add_argument('-g', '--access',
                    help="""Sequencing-accessible genomic regions, or exons to
                    use as possible targets (e.g. output of refFlat2bed.py)""")
    AP.add_argument('-t', '--targets',
                    help="""Potentially targeted genomic regions, e.g. all
                    possible exons for the reference genome. Format: BED,
                    interval list, etc.""")
    AP.add_argument('-o', '--output', help="Output filename.")
    main(AP.parse_args())
