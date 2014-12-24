#!/usr/bin/env python

"""Convert THetA output to a BED-like, CNVkit-like tabular format.

Equivalently, use the THetA results file to convert CNVkit .cns segments to
integer copy number calls.
"""
from __future__ import division, print_function

import math

import cnvlib

PLOIDY = 2

def parse_theta_results(fname):
    """Parse THetA results into a data structure.

    Columns: NLL, mu, C, p*
    """
    with open(fname) as handle:
        header = next(handle).rstrip().split('\t')
        body = next(handle).rstrip().split('\t')
        assert len(body) == len(header) == 4

        # NLL
        nll = float(body[0])

        # mu
        mu = body[1].split(',')
        mu_normal = float(mu[0])
        mu_tumors = map(float, mu[1:])

        # C
        copies = body[2].split(':')
        if len(mu_tumors) == 1:
            # 1D array of integers
            # Replace X with None for "missing"
            copies = [[int(c) if c.isdigit() else None
                       for c in copies]]
        else:
            # List of lists of integer-or-None (usu. 2 x #segments)
            copies = [[int(c) if c.isdigit() else None
                       for c in subcop]
                      for subcop in zip(*[c.split(',') for c in copies])]

        # p*
        probs = body[3].split(',')
        if len(mu_tumors) == 1:
            # 1D array of floats, or None for "X" (missing/unknown)
            probs = [float(p) if not p.isalpha() else None
                     for p in probs]
        else:
            probs = [[float(p) if not p.isalpha() else None
                      for p in subprob]
                     for subprob in zip(*[p.split(',') for p in probs])]
    return {"NLL": nll,
            "mu_normal": mu_normal,
            "mu_tumors": mu_tumors,
            "C": copies,
            "p*": probs}


def main(args):
    """."""
    tumor_segs = cnvlib.read(args.tumor_cns)
    theta = parse_theta_results(args.theta_results)
    for i, copies in enumerate(theta['C']):
        # Replace segment values with these integers
        # Drop any segments where the C value is None
        new_segs = []
        for seg, ncop in zip(tumor_segs.copy(), copies):
            if ncop is None:
                continue
            seg["coverage"] = math.log((ncop or 0.5) / PLOIDY, 2)
            new_segs.append(seg)
        new_cns = tumor_segs.to_rows(new_segs)
        new_cns.write("%s-%d.cni" % (tumor_segs.sample_id, i + 1))


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("tumor_cns")
    AP.add_argument("theta_results")
    main(AP.parse_args())
