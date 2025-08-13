#!/usr/bin/env python
"""Extract target and antitarget BED files from a CNVkit reference file.

Once you have a stable CNVkit reference for your platform, you can use this
script to drop the "bad" bins from your target and antitarget BED files and
avoid unnecessarily calculating coverage in those bins during future runs.

This script is also useful to recover the target and antitarget BED files that
match the reference if those BED files are missing or you're not sure which ones
are correct.
"""
import argparse
import logging

from .. import reference, read
from skgenome import tabio


def argument_parsing():
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("reference", help="Reference file.")
    AP.add_argument("-o", "--output",
                    help="Output base name (extensions added automatically).")
    return AP.parse_args()


def reference2targets(args) -> None:
    ref = read(args.reference)
    targets, antitargets = reference.reference2regions(ref)
    name = args.output or ref.sample_id
    tabio.write(targets, name + '.target.bed', 'bed4')
    tabio.write(antitargets, name + '.antitarget.bed', 'bed4')


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    arguments = argument_parsing()
    reference2targets(arguments)


if __name__ == '__main__':
    main()
