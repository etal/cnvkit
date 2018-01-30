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
                 "from "+ args.in_fmt if args.in_fmt != 'auto' else '',
                 args.out_fmt)

    # TODO - add back merge/flatten/exon options from refFlat2bed
    table = tabio.read(args.infile, args.in_fmt)
    tabio.write(table, args.output, args.out_fmt)



if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('infile',
                    default=sys.stdin, nargs='?', metavar="FILE",
                    help="Input filename. [Default: stdin]")
    AP.add_argument('-f', '--from',
                    default='auto', dest='in_fmt', metavar="FORMAT",
                    help="Input format. [Default: auto-detect]")
    AP.add_argument('-t', '--to',
                    required=True, dest='out_fmt', metavar="FORMAT",
                    help="Output format. [Required]")
    AP.add_argument('-o', '--output', metavar="FILE",
                    help="Output filename. [Default: stdout]")
    main(AP.parse_args())
