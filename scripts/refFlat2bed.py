#!/usr/bin/env python

"""Extract canonical gene coordinates from UCSC's refFlat.txt.

Usage:
    refFlat2bed.py /path/to/refFlat.txt > genes.bed
    refFlat2bed.py -e /path/to/refFlat.txt > exons.bed
"""
from __future__ import print_function

import argparse
import sys

from cnvlib import tabio


AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument('refflat',
                help="UCSC refFlat.txt for the reference genome.")
AP.add_argument('-e', '--exons', action='store_true',
                help="""Emit each exon instead of the whole gene regions.""")
AP.add_argument('-f', '--flatten', action='store_true',
                help="""Flatten overlapping regions, keeping original
                boundaries. Not recommended with --exons.""")
AP.add_argument('-m', '--merge',
                metavar='BASEPAIRS', nargs='?', type=int, const=1,
                help="""Merge overlapping regions with different names.
                Recommended with --exons. Optional argument value is the
                number of overlapping bases between two regions to trigger a
                merge. [Default: %(const)s]""")
AP.add_argument('-o', '--output',
                help="Output filename.")
args = AP.parse_args()

regions = tabio.read(args.refflat, 'refflat', exons=args.exons)
if args.flatten:
    regions = regions.flatten()
elif args.merge:
    regions = regions.merge(bp=args.merge)
tabio.write(regions, args.output or sys.stdout, 'bed4')
