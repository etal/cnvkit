"""Import from other formats to the CNVkit format."""
from __future__ import absolute_import, division, print_function
from builtins import map, next, zip

import logging

import numpy as np
from skgenome import tabio

from . import params


# __________________________________________________________________________
# import-picard

def do_import_picard(fname, too_many_no_coverage=100):
    garr = tabio.read(fname, "picardhs")
    garr["gene"] = garr["gene"].apply(unpipe_name)
    # Create log2 column from coverages, avoiding math domain error
    coverages = garr["ratio"].copy()
    no_cvg_idx = (coverages == 0)
    if no_cvg_idx.sum() > too_many_no_coverage:
        logging.warning("WARNING: Sample %s has >%d bins with no coverage",
                        garr.sample_id, too_many_no_coverage)
    coverages[no_cvg_idx] = 2**params.NULL_LOG2_COVERAGE
    garr["log2"] = np.log2(coverages)
    return garr


def unpipe_name(name):
    """Fix the duplicated gene names Picard spits out.

    Return a string containing the single gene name, sans duplications and pipe
    characters.

    Picard CalculateHsMetrics combines the labels of overlapping intervals
    by joining all labels with '|', e.g. 'BRAF|BRAF' -- no two distinct
    targeted genes actually overlap, though, so these dupes are redundant.
    Meaningless target names are dropped, e.g. 'CGH|FOO|-' resolves as 'FOO'.
    In case of ambiguity, the longest name is taken, e.g. "TERT|TERT Promoter"
    resolves as "TERT Promoter".
    """
    if '|' not in name:
        return name
    gene_names = set(name.split('|'))
    if len(gene_names) == 1:
        return gene_names.pop()
    cleaned_names = gene_names.difference(params.IGNORE_GENE_NAMES)
    if cleaned_names:
        gene_names = cleaned_names
    new_name = sorted(gene_names, key=len, reverse=True)[0]
    if len(gene_names) > 1:
        logging.warning("WARNING: Ambiguous gene name %r; using %r",
                        name, new_name)
    return new_name


# __________________________________________________________________________
# import-theta

def do_import_theta(segarr, theta_results_fname, ploidy=2):
    theta = parse_theta_results(theta_results_fname)
    # THetA doesn't handle sex chromosomes well
    segarr = segarr.autosomes()
    for copies in theta['C']:
        if len(copies) != len(segarr):
            copies = copies[:len(segarr)]
        # Drop any segments where the C value is None
        mask_drop = np.array([c is None for c in copies], dtype='bool')
        segarr = segarr[~mask_drop].copy()
        ok_copies = np.asfarray([c for c in copies if c is not None])
        # Replace remaining segment values with these integers
        segarr["cn"] = ok_copies.astype('int')
        ok_copies[ok_copies == 0] = 0.5
        segarr["log2"] = np.log2(ok_copies / ploidy)
        segarr.sort_columns()
        yield segarr


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
        mu_tumors = list(map(float, mu[1:]))

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
