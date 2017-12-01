#!/usr/bin/env python
"""Convert a cohort of per-gene log2 ratios to CNVkit .cnr format.
"""
from __future__ import absolute_import, division, print_function

import logging
import os
import sys

import pandas as pd
from skgenome import tabio

from cnvlib import rna


def aggregate_gene_counts(filenames):
    prev_row_count = None
    sample_cols = {}
    for fname in filenames:
        d = (pd.read_table(fname,
                           header=None,
                           comment="_",
                           names=["gene_id", "expected_count"],
                           converters={"gene_id": rna.before(".")})
             .set_index("gene_id"))
             # .drop(["__no_feature", "__ambiguous", "__too_low_aQual",
             # "__not_aligned", "__alignment_not_unique"]))
        if prev_row_count is None:
            prev_row_count = len(d)
        elif len(d) != prev_row_count:
            raise RuntimeError("Number of rows in each input file is not equal")
        sample_id = rna.before(".")(os.path.basename(fname))
        sample_cols[sample_id] = d.expected_count.fillna(0)
    sample_counts = pd.DataFrame(sample_cols)
    return sample_counts


def main(args):
    sample_counts = aggregate_gene_counts(args.gene_counts)
    sample_counts = rna.filter_probes(sample_counts)
    # DBG
    if args.output:
        sample_counts.to_csv(args.output + ".sample_counts.tsv",
                             sep='\t', index=True)
        print("Wrote", args.output + ".sample_counts.tsv",
              "with", len(sample_counts), "rows")

    if args.tcga_cnv and args.tcga_expression:
        logging.info("Loading gene metadata "
                     "and TCGA gene expression/CNV profiles")
    else:
        logging.info("Loading gene metadata")
    gene_info = rna.load_gene_info(args.gene_resource,
                                   args.tcga_cnv, args.tcga_expression)

    print("Aligning gene info to sample gene counts")
    gene_info, sample_counts, sample_data_log2 = rna.align_gene_info_to_samples(
        gene_info, sample_counts, None)

    print("Writing output files")
    # Summary table has log2-normalized values, not raw counts
    # ENH show both, with column header suffixes to distinguish?
    all_data = pd.concat([gene_info, sample_data_log2], axis=1)
    if args.output:
        all_data.to_csv(args.output, sep='\t', index=True)
        print("Wrote", args.output, "with", len(all_data), "rows")
    else:
        print(all_data.describe(), file=sys.stderr)

    if args.cnr_dir:
        # CNVkit files have both absolute and log2-normalized read counts
        cnrs = rna.attach_gene_info_to_cnr(sample_counts, sample_data_log2,
                                          gene_info)
        for cnr in cnrs:
            outfname = os.path.join(args.cnr_dir, cnr.sample_id + ".cnr")
            cnr = rna.correct_cnr(cnr)
            tabio.write(cnr, outfname, 'tab')



if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("gene_counts", nargs='+',
                    help="""Tabular files with 2 columns: Ensembl gene ID and
                    number of reads mapped to each gene.""")
    AP.add_argument("-o", "--output", metavar="FILE",
                    help="Output file name (summary table).")
    AP.add_argument("--cnr-dir", metavar="PATH",
                    help="""Directory to write a CNVkit .cnr file for each input
                    sample.""")
    AP.add_argument("-g", "--gene-resource", metavar="FILE",
                    help="Location of gene info table from BioMart.")

    AP_tcga = AP.add_argument_group(
        "To correlate each gene's copy number with expression")
    AP_tcga.add_argument("--tcga-cnv", metavar="CNV_FILE",
                    help="""CNVs for many TCGA samples, downloaded from
                    cBioPortal.""")
    AP_tcga.add_argument("--tcga-expression", metavar="RNA_FILE",
                    help="""Gene expression for many TCGA samples, downloaded from
                    cBioPortal.""")

    main(AP.parse_args())
