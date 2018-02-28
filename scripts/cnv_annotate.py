#!/usr/bin/env python

"""Update gene names in CNVkit .cnn/.cnr files.
"""
from __future__ import absolute_import, division, print_function
import argparse
import logging
import sys

from skgenome import tabio
from cnvlib.cmdutil import read_cna

logging.basicConfig(level=logging.INFO, format="%(message)s")


def main(args):
    annot = tabio.read_auto(args.annotate)
    cnarr = read_cna(args.cn_file)
    cnarr['gene'] = annot.into_ranges(cnarr, 'gene', '-')
    tabio.write(cnarr, args.output or sys.stdout)
    # ENH:
    #  --short-names
    #  --no-antitarget
    #  annotation: --merge, --flatten, --exons, ...
    #  cut antitargets & insert untargeted gene names
    #      some math for how to update probes, weight



if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('annotate', help="Genome annotations.")
    AP.add_argument('cnv_file', help="CNVkit .cnn or .cnr file.")
    AP.add_argument('-o', '--output', help="Output filename.")
    main(AP.parse_args())
