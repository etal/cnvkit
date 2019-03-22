#!/usr/bin/env python
"""Z-test for single-bin copy number alterations."""
import argparse
import logging
import sys

import cnvlib
from cnvlib.cmdutil import verify_sample_sex
from cnvlib.ztest import do_ztest
from skgenome import tabio

logging.basicConfig(level=logging.INFO, format="%(message)s")


def _cmd_ztest(args):
    cnarr = cnvlib.read(args.cnarr)
    if args.segment:
        segments = cnvlib.read(args.segment)
        is_sample_female = None
    else:
        segments = None
        is_sample_female = verify_sample_sex(cnarr, args.sample_sex,
                                             args.male_reference)
    sig = do_ztest(cnarr, segments, args.male_reference, is_sample_female,
                   args.alpha, args.target)
    if len(sig):
        tabio.write(sig, args.output or sys.stdout)



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("cnarr", help=".cnr file")
    AP.add_argument('-s', '--segment', metavar="FILENAME",
                    help="""Segmentation calls (.cns), the output of the
                    'segment' command).""")

    AP.add_argument("-a", "--alpha", type=float, default=0.005,
                    help="Significance threhold. [Default: %(default)s]")
    AP.add_argument("-t", "--target", action="store_true",
                    help="Test target bins only; ignore off-target bins.")
    AP.add_argument('-y', '--male-reference', '--haploid-x-reference',
                    action='store_true',
                    help="""Assume inputs were normalized to a male reference
                    (i.e. female samples will have +1 log-coverage of chrX;
                    otherwise male samples would have -1 chrX).""")
    AP.add_argument('-x', '--sample-sex',
                    choices=('m', 'y', 'male', 'Male',
                             'f', 'x', 'female', 'Female'),
                    help="""Specify the sample's chromosomal sex as male or
                    female.  (Otherwise guessed from X and Y coverage).""")

    AP.add_argument("-o", "--output",
                    help="Output filename.")
    _cmd_ztest(AP.parse_args())
