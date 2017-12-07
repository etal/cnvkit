#!/usr/bin/env python
"""Extract CNVkit-compatible per-gene, per-sample log2 ratios from RSEM outputs.
"""
from __future__ import absolute_import, division, print_function

import logging
import os
import sys

import numpy as np
import pandas as pd
from skgenome import tabio

from cnvlib import rna


def aggregate_rsem(fnames):
    """Pull out the expected read counts from each RSEM file.

    The expected read counts are located in the fifth column of this version of
    RSEM output. They start on the second line (after a header line). I will
    also pull the ensembl gene id, which is located in the first column.

    Returns:
        sample_counts : DataFrame
            Row index is Ensembl gene ID, column index is filename.
        tx_lengths : Series
            Gene lengths.
    """
    prev_row_count = None
    sample_cols = {}
    length_cols = []
    length_colname = 'length'  # or: 'effective_length'
    for fname in fnames:
        # NB: read_table(index_col=_) works independently of combine=, dtype=
        #   so index column needs to be processed separately
        #   https://github.com/pandas-dev/pandas/issues/9435
        d = pd.read_table(fname,
                          usecols=['gene_id', length_colname, 'expected_count'],
                          #  index_col='gene_id',
                          converters={'gene_id': rna.before('.')}
                         ).set_index('gene_id')
        if prev_row_count is None:
            prev_row_count = len(d)
        elif len(d) != prev_row_count:
            raise RuntimeError("Number of rows in each input file is not equal")
        sample_id = rna.before(".")(os.path.basename(fname))
        sample_cols[sample_id] = d.expected_count.fillna(0)
        length_cols.append(d[length_colname])
    sample_counts = pd.DataFrame(sample_cols)
    tx_lengths = pd.Series(np.vstack(length_cols).mean(axis=0),
                             index=sample_counts.index)
    return sample_counts, tx_lengths


def main(args):
    print("Loading", len(args.rsem_results), "RSEM files")
    sample_counts, tx_lengths = aggregate_rsem(args.rsem_results)

    sample_counts = rna.filter_probes(sample_counts)
    # DBG
    if args.output:
        sample_counts.to_csv(args.output + ".sample_counts.tsv",
                             sep='\t', index=True)
        print("Wrote", args.output + ".sample_counts.tsv",
              "with", len(sample_counts), "rows")

    print("Loading gene metadata",
          "and TCGA gene expression/CNV profiles" if args.correlations else "")
    gene_info = rna.load_gene_info(args.gene_resource, args.correlations)

    print("Aligning gene info to sample gene counts")
    gene_info, sample_counts, sample_depths_log2 = rna.align_gene_info_to_samples(
        gene_info, sample_counts, tx_lengths)

    print("Writing output files")
    # Summary table has log2-normalized values, not raw counts
    # ENH show both, with column header suffixes to distinguish?
    all_data = pd.concat([gene_info, sample_depths_log2], axis=1)
    if args.output:
        all_data.to_csv(args.output, sep='\t', index=True)
        print("Wrote", args.output, "with", len(all_data), "rows")
    else:
        print(all_data.describe(), file=sys.stderr)

    if args.cnr_dir:
        # CNVkit files have both absolute and log2-normalized read counts
        cnrs = rna.attach_gene_info_to_cnr(sample_counts, sample_depths_log2,
                                           gene_info)
        for cnr in cnrs:
            outfname = os.path.join(args.cnr_dir, cnr.sample_id + ".cnr")
            cnr = rna.correct_cnr(cnr)
            tabio.write(cnr, outfname, 'tab')



if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("rsem_results", nargs='+',
                    help="RSEM output files (*.results)")
    AP.add_argument("-o", "--output", metavar="FILE",
                    help="Output file name (summary table).")
    AP.add_argument("-d", "--cnr-dir", metavar="PATH",
                    help="""Directory to write a CNVkit .cnr file for each input
                    sample.""")
    AP.add_argument("-g", "--gene-resource", metavar="FILE",
                    help="Location of gene info table from BioMart.")
    AP.add_argument("-c", "--correlations", metavar="FILE",
                    help="""Correlation of each gene's copy number with
                    expression. Output of cnv_expression_correlate.py.""")

    main(AP.parse_args())
