#!/usr/bin/env python
"""[DEPRECATED; see 'autobin'] Calculate good average bin sizes from a BAM."""
from __future__ import absolute_import, division, print_function

import argparse
import logging

from cnvlib import autobin, tabio


logging.basicConfig(level=logging.WARN, format="%(message)s")

row_labels = {
    "amplicon": ("Targets (sampling)", ""),
    "hybrid": ("On-target", "Off-target"),
    "wgs": ("Genome", ""),
}


def read_regions(bed_fname):
    if bed_fname:
        regions = tabio.read_auto(bed_fname)
        if len(regions):
            return regions
        else:
            logging.warn("No regions to estimate depth from %s",
                         regions.meta.get('filename', ''))
    return None


def print_table(fields, labels):
    """Print label-depth pairs as a table, with calculated bin size."""
    width = max(map(len, labels)) + 1
    print(' ' * width, "Depth", "Bin size", sep='\t')
    for label, (depth, binsize) in zip(labels, fields):
        if depth is not None:
            print((label + ":").ljust(width),
                  format(depth, '.3f'),
                  binsize,
                  sep='\t')


def main(args):
    if args.method in ('hybrid', 'amplicon') and not args.targets:
        raise RuntimeError("Sequencing method %r requires targets", args.method)
    elif args.method == 'wgs' and args.targets:
        logging.warn("Targets will be ignored: %s", args.targets)
    if args.method == 'amplicon' and args.access:
        logging.warn("Sequencing-accessible regions will be ignored: %s",
                        args.access)

    access = read_regions(args.access)
    targets = read_regions(args.targets)
    fields = autobin.do_autobin(args.bam, args.method, targets, access, args.bp_per_bin)
    print_table(fields, labels=row_labels[args.method])


if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('bam',
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
    AP.add_argument('-b', '--bp-per-bin', type=float, default=100000.,
                    help="""Desired average number of bases per bin.
                    [Default: %(default)s]""")
    AP.add_argument('-o', '--output', help="Output filename.")
    main(AP.parse_args())
