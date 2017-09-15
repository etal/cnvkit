#!/usr/bin/env python
"""Z-test for single-bin copy number alterations."""
from __future__ import division, print_function

import argparse
import logging
import sys

import numpy as np
from scipy.stats import norm

import cnvlib
from skgenome import tabio


def ztest(cnarr, alpha, verbose=True):
    """Get a probability for each bin based on its Z-score.

    Return bins where the probability < `alpha`.
    """
    # Bin weights ~ 1/variance; bin log2 values already centered at 0.0
    sd = np.sqrt(1 - cnarr['weight'])
    p = norm.pdf(cnarr['log2'], loc=0, scale=sd)
    # Correct for multiple hypothesis tests
    # p *= len(p)       # Bonferroni
    p = p_adjust_bh(p)  # B-H
    cnarr['ztest'] = p

    is_sig = cnarr['ztest'] < alpha
    if verbose:
        print("Significant hits in", is_sig.sum(), "/", len(is_sig), "bins",
              "({:.4g}%)".format(100 * is_sig.sum() / len(is_sig)))
    return cnarr[is_sig]


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
                    help="Significance threhold")
    AP.add_argument("-t", "--target", action="store_true",
                    help="Test target bins only; ignore off-target bins.")
    AP.add_argument("-o", "--output",
                    help="Output filename.")
    args = AP.parse_args()

    cnarr = cnvlib.read(args.cnarr)

    if args.target:
        antitarget_idx = cnarr['gene'].isin(cnvlib.params.ANTITARGET_ALIASES)
        if antitarget_idx.any():
            print("Ignoring", antitarget_idx.sum(), "off-target bins")
            cnarr = cnarr[~antitarget_idx]

    sig = ztest(cnarr, args.alpha)
    if len(sig):
        tabio.write(sig, args.output or sys.stdout)
