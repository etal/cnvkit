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

    table = tabio.read(args.infile, args.in_fmt)
    tabio.write(table, args.output, args.out_fmt)



if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('infile', nargs='?', default=sys.stdin,
                    help="Input filename.")
    AP.add_argument('-f', '--from', dest='in_fmt', default='auto',
                    help="Input format.")
    AP.add_argument('-t', '--to', dest='out_fmt', required=True,
                    help="Output format.")
    AP.add_argument('-o', '--output',
                    help="Output filename.")
    main(AP.parse_args())
