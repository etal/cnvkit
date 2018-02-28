#!/usr/bin/env python
"""Z-test for single-bin copy number alterations."""
from __future__ import division, print_function

import argparse
import logging
import sys

import numpy as np
from scipy.stats import norm

import cnvlib
from cnvlib.cmdutil import verify_sample_sex
from skgenome import tabio


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


def do_ztest(cnarr, segments=None, is_male_reference=False,
             is_sample_female=None, alpha=0.005, target_only=False):
    """Get a probability for each bin based on its Z-score.

    Return bins where the probability < `alpha`.
    """
    if segments:
        # Subtract segment means to report only the CNA bins that weren't
        # already detected (including exon-size CNAs within a larger-scale,
        # smaller-amplitude CNA)
        cnarr['log2'] = cnarr.residuals(segments)
    else:
        # Account for non-diploid sex chromosomes
        cnarr['log2'] -= cnarr.expect_flat_log2(is_male_reference)
        if is_sample_female:
            # chrX has same ploidy as autosomes; chrY is just unusable noise
            cnarr = cnarr[cnarr.chromosome != cnarr._chr_y_label]
        else:
            # 1/2 #copies of each sex chromosome
            is_chr_xy = cnarr.chromosome.isin((cnarr._chr_x_label,
                                               cnarr._chr_y_label))
            cnarr[is_chr_xy, 'log2'] += 1.0

    if target_only:
        antitarget_idx = cnarr['gene'].isin(cnvlib.params.ANTITARGET_ALIASES)
        if antitarget_idx.any():
            print("Ignoring", antitarget_idx.sum(), "off-target bins")
            cnarr = cnarr[~antitarget_idx]

    cnarr['ztest'] = z_prob(cnarr)
    is_sig = cnarr['ztest'] < alpha
    logging.info("Significant hits in {}/{} bins ({:.3g}%)"
                 .format(is_sig.sum(), len(is_sig),
                         100 * is_sig.sum() / len(is_sig)))
    return cnarr[is_sig]


def z_prob(cnarr):
    # Bin weights ~ 1-variance; bin log2 values already centered at 0.0
    sd = np.sqrt(1 - cnarr['weight'])
    # Convert to Z-scores
    z = cnarr['log2'] / sd
    # Two-sided survival function (1-CDF) probability
    p = 2. * norm.cdf(z)
    # Similar to the above -- which is better?
    # p2 = 2 * norm.pdf(cnarr['log2'], loc=0, scale=sd)
    # if not np.allclose(p, p2):
    #     print("Max diff:", np.abs(p - p2).max())
    #     print("Median diff:", np.median(np.abs(p - p2)))
    #     print("Ratio:", (p / p2).mean())
    # Correct for multiple hypothesis tests
    return p_adjust_bh(p)


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("cnarr", help=".cnr file")
    AP.add_argument("-a", "--alpha", type=float, default=0.005,
                    help="Significance threhold. [Default: %(default)s]")
    AP.add_argument("-t", "--target", action="store_true",
                    help="Test target bins only; ignore off-target bins.")

    AP.add_argument('-s', '--segment', metavar="FILENAME",
                    help="""Segmentation calls (.cns), the output of the
                    'segment' command).""")
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
