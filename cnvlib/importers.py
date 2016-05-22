"""Import from other formats to the CNVkit format."""
from __future__ import absolute_import, division, print_function
from builtins import map, next, zip

import os.path
import subprocess


# __________________________________________________________________________
# import-picard

def find_picard_files(file_and_dir_names):
    """Search the given paths for 'targetcoverage' CSV files.

    Per the convention we use in our Picard applets, the target coverage file
    names end with '.targetcoverage.csv'; anti-target coverages end with
    '.antitargetcoverage.csv'.
    """
    filenames = []
    for tgt in file_and_dir_names:
        if os.path.isdir(tgt):
            # Collect the target coverage files from this directory tree
            fnames = subprocess.check_output(['find', tgt,
                                              '-name', '*targetcoverage.csv']
                                            ).splitlines()
            if not fnames:
                raise RuntimeError("Given directory %s does not contain any "
                                   "'*targetcoverage.csv' files."
                                   % tgt)
            filenames.extend(fnames)
        elif os.path.isfile(tgt):
            filenames.append(tgt)
        else:
            raise ValueError("Given path is neither a file nor a directory: %s"
                             % tgt)
    filenames.sort()
    return filenames


# __________________________________________________________________________
# import-theta

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
