#!/usr/bin/env python

"""Extract target and antitarget BED files from a CNVkit reference file.

Once you have a stable CNVkit reference for your platform, you can use this
script to drop the "bad" bins from your target and antitarget BED files and
avoid unnecessarily calculating coverage in those bins during future runs.

This script is also useful to recover the target and antitarget BED files that
match the reference if those BED files are missing or you're not sure which ones
are correct.
"""

import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")

import cnvlib
from cnvlib import ngfrills, reference


def write_bed(rows, fname):
    """Write region coordinates to `fname` in BED format."""
    with ngfrills.safe_write(fname, False) as outfile:
        i = 0
        for i, row in enumerate(rows):
            outfile.write("\t".join(map(str, row)) + '\n')
        logging.info("Wrote %s with %d bins", fname, i + 1)


def main(args):
    """Run the script."""
    ref = cnvlib.read(args.reference)
    targets, antitargets = reference.reference2regions(ref)
    name = args.output or ref.sample_id
    write_bed(targets, name + '.target.bed')
    write_bed(antitargets, name + '.antitarget.bed')


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("reference", help="Reference file.")
    AP.add_argument("-o", "--output",
                    help="Output base name (extensions added automatically).")
    main(AP.parse_args())

