#!/usr/bin/env python
"""Update gene names in CNVkit .cnn/.cnr files.
"""
import argparse
import logging
import sys

from skgenome import tabio
from ..cmdutil import read_cna


def argument_parsing():
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('annotate', help="Genome annotations.")
    AP.add_argument('cnv_file', help="CNVkit .cnn or .cnr file.")
    AP.add_argument('-o', '--output', help="Output filename.")
    return  AP.parse_args()


def cnv_annotate(args) -> None:
    annot = tabio.read_auto(args.annotate)
    cnarr = read_cna(args.cnv_file)
    cnarr['gene'] = annot.into_ranges(cnarr, 'gene', '-')
    tabio.write(cnarr, args.output or sys.stdout)
    # ENH:
    #  --short-names
    #  --no-antitarget
    #  annotation: --merge, --flatten, --exons, ...
    #  cut antitargets & insert untargeted gene names
    #      some math for how to update probes, weight


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    arguments = argument_parsing()
    cnv_annotate(arguments)


if __name__ == '__main__':
    main()
