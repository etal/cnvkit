#!/usr/bin/env python
"""Convert between tabular formats using scikit-genome I/O."""
from __future__ import absolute_import, division, print_function
import argparse
import logging
import sys

from skgenome import tabio

logging.basicConfig(level=logging.INFO, format="%(message)s")


def main(args):
    logging.info("Converting %s%s to %s",
                 "input" if args.infile is sys.stdin else args.infile,
                 " from "+ args.in_fmt if args.in_fmt != 'auto' else '',
                 args.out_fmt)

    if args.in_fmt == 'auto':
        args.in_fmt = tabio.sniff_region_format(args.infile)
    # Format-specific input options
    kwargs = {}
    if args.in_fmt == 'gff':
        if args.gff_tag:
            kwargs['tag'] = args.gff_tag
        if args.gff_type:
            kwargs['keep_type'] = args.gff_type
    elif args.in_fmt == 'refflat':
        if args.refflat_type == 'exon':
            kwargs['exons'] = True
        elif args.refflat_type == 'cds':
            kwargs['cds'] = True
    regions = tabio.read(args.infile, args.in_fmt, **kwargs)

    # Post-processing
    if args.flatten:
        regions = regions.flatten()
    elif args.merge:
        regions = regions.merge(bp=args.merge)

    tabio.write(regions, args.output, args.out_fmt)


if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('infile', metavar="FILE",
                    default=sys.stdin, nargs='?',
                    help="Input filename. [Default: stdin]")
    AP.add_argument('-o', '--output', metavar="FILE",
                    help="Output filename. [Default: stdout]")
    AP.add_argument('-f', '--from', metavar="FORMAT",
                    default='auto', dest='in_fmt',
                    help="Input format. [Default: auto-detect]")
    AP.add_argument('-t', '--to', metavar="FORMAT",
                    default='bed4', dest='out_fmt',
                    help="Output format. [Required]")

    AP.add_argument('--flatten', action='store_true',
                    help="""Flatten overlapping regions, keeping original
                    boundaries. Not recommended with --exons.""")
    AP.add_argument('--merge', metavar='BASEPAIRS',
                    nargs='?', type=int, const=1,
                    help="""Merge overlapping regions with different names.
                    Recommended with --exons. Optional argument value is the
                    number of overlapping bases between two regions to trigger a
                    merge. [Default: %(const)s]""")
    # ENH combine --gff-type, --refflat-type
    # AP.add_argument('-e', '--exons', action='store_true',
    #                 help="""Emit each exon instead of the whole gene regions.
    #                 (Applies to 'gff' and 'refflat' input formats.)""")

    AP_fmt = AP.add_argument_group("Format-specific options")
    AP_fmt.add_argument("--gff-tag",
                    help="""GFF attributes tag to use for gene names.""")
    AP_fmt.add_argument("--gff-type",
                    help="""GFF record type (column 3) to use exclusively.
                    E.g. for GFF3 files, 'exon', 'gene', and other SOFA terms
                    can be used here.""")
    AP_fmt.add_argument('--refflat-type',
                    choices=('exon', 'cds', 'transcript'), default='transcript',
                    help="""Emit each exon instead of the whole gene regions.""")

    main(AP.parse_args())
