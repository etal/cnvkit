"""Import from other formats to the CNVkit format."""
from __future__ import absolute_import, division, print_function
from builtins import next
from builtins import map
from builtins import zip

import math
import os.path
import subprocess

import pandas as pd

from .cnary import CopyNumArray as CNA


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
# import-seg

LOG2_10 = math.log(10, 2)   # To convert log10 values to log2

def import_seg(segfname, chrom_names, chrom_prefix, from_log10):
    """Parse a SEG file as an iterable of CopyNumArray instances.

    `chrom_names`:
        Map (string) chromosome IDs to names. (Applied before chrom_prefix.)
        e.g. {'23': 'X', '24': 'Y', '25': 'M'}

    `chrom_prefix`: prepend this string to chromosome names
        (usually 'chr' or None)

    `from_log10`: Convert values from log10 to log2.
    """
    dframe = pd.read_table(segfname, na_filter=False)
    if len(dframe.columns) == 6:
        dframe.columns = ['sample_id', 'chromosome', 'start', 'end', 'nprobes',
                          'mean']
    elif len(dframe.columns) == 5:
        dframe.columns = ['sample_id', 'chromosome', 'start', 'end', 'mean']
    else:
        raise ValueError("SEG format expects 5 or 6 columns; found {}: {}"
                         .format(len(dframe.columns), ' '.join(dframe.columns)))

    # Calculate values for output columns
    dframe['chromosome'] = dframe['chromosome'].apply(str)
    if chrom_names:
        dframe['chromosome'] = dframe['chromosome'].apply(lambda c:
                                                          chrom_names.get(c, c))
    if chrom_prefix:
        dframe['chromosome'] = dframe['chromosome'].apply(lambda c:
                                                          chrom_prefix + c)
    if from_log10:
        dframe['mean'] *= LOG2_10
    dframe['gene'] = ["G" if mean >= 0 else "L" for mean in dframe['mean']]

    for sid in pd.unique(dframe['sample_id']):
        sample = dframe[dframe['sample_id'] == sid]
        cols = {'chromosome': sample['chromosome'],
                'start': sample['start'],
                'end': sample['end'],
                'gene': sample['gene'],
                'log2': sample['mean']}
        if 'nprobes' in dframe:
            cols['probes'] = sample['nprobes']
        cns = CNA.from_columns(cols, {'sample_id': sid})
        cns.sort()
        yield cns


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
