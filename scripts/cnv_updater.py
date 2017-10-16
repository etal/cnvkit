#!/usr/bin/env python
"""Update .cnn/.cnr files from older CNVkit versions to match current defaults.

CNVkit v0.8.0 and later uses a 'depth' column in the *.targetcoverage.cnn and
*.antitargetcoverage.cnn files produced by the 'coverage' command.
To use .cnn files created by CNVkit v0.7.11 or earlier with the current version,
run this script on the old .cnn files to convert them to the new format,
calculating 'depth' from 'log2'.

CNVkit v0.9.0 and later automatically names off-target bins "Antitarget",
instead of the previous "Background". This script also updates those bin names
in the 'gene' column.
"""
from __future__ import absolute_import, division, print_function

import argparse
import logging
import os.path

import numpy as np

import cnvlib
from skgenome import tabio

logging.basicConfig(level=logging.INFO, format="%(message)s")


AP = argparse.ArgumentParser(description=__doc__)
AP.add_argument("cnn_files", nargs='+',
        help="""CNVkit coverage files to update (*.targetcoverage.cnn,
                *.antitargetcoverage.cnn).""")
AP.add_argument("-d", "--output-dir", default=".",
        help="""Directory to write output .cnn files.""")
AP.add_argument("-s", "--suffix", default=".updated",
        help="""Filename suffix to add before the '.cnn' extension in output
                files. [Default: %(default)s]""")
args = AP.parse_args()

for fname in args.cnn_files:
    cnarr = cnvlib.read(fname)
    # Convert coverage depths from log2 scale to absolute scale.
    # NB: The log2 values are un-centered in CNVkit v0.7.0(?) through v0.7.11;
    # earlier than that, the average 'depth' will be about 1.0.
    cnarr['depth'] = np.exp2(cnarr['log2'])
    # Rename "Background" bins to "Antitarget"
    # NB: The default off-target bin name was changed in CNVkit v0.9.0
    cnarr['gene'] = cnarr['gene'].replace("Background",
                                          cnvlib.params.ANTITARGET_NAME)
    cnarr.sort_columns()
    # Construct the output filename
    base, ext = os.path.basename(fname).rsplit('.', 1)
    if '.' in base:
        base, zone = base.rsplit('.', 1)
        out_fname = '.'.join((base + args.suffix, zone, ext))
    else:
        # e.g. reference.cnn or .cnr file, no "*.targetcoverage.*" in name
        out_fname = '.'.join((base + args.suffix, ext))
    tabio.write(cnarr, os.path.join(args.output_dir, out_fname))
