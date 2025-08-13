#!/usr/bin/env python
"""Calculate correlation coefficients for gene expression and copy number.

Data source for both inputs is TCGA via cBioPortal.
"""
import logging
import sys
import warnings
import argparse

import pandas as pd
from scipy.stats import spearmanr, kendalltau

from ..rna import before

def argument_parsing():
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument("cnv_fname", metavar="CNV_FILE",
                    help="""Gene copy number calls for many samples.""")
    AP.add_argument("expression_fname", metavar="RNA_FILE",
                    help="""Gene expression for many samples (mostly
                    overlapping with CNV samples).""")
    AP.add_argument("-o", "--output", metavar="FILE",
                    help="Output file name (summary table).")
    return AP.parse_args()

def correlate_cnv_expression(cnv_fname, expression_fname):
    """Get correlation coefficients for matched copy number and expression data.

    cBioPortal offers a nice feature in which you can download a summary of many
    large-scale sequencing studies. In this summary are two files that contain
    the copy number and expression values of every gene in the study for every
    sample.  This summary is available for nearly every TCGA study, and the data
    is intuitive to access, therefore I have designed this pre-processing script
    to accept these as inputs. Of course, the user can calculate their own
    Pearson values from other sources of data if they prefer -- in this case,
    the user should formate their data to match the output of this prepocessing
    script.
    """
    shared_key = 'Entrez_Gene_Id'
    cnv_table = load_tcga_table(cnv_fname, shared_key)
    expr_table = load_tcga_table(expression_fname, shared_key)
    # Ensure rows and columns match
    cnv_table, expr_table = cnv_table.align(expr_table, join='inner')
    logging.info("Trimmed TCGA tables to shape: %s", cnv_table.shape)

    # Calculate correlation coefficients
    logging.info("Calculating correlation coefficients")
    c_nums = cnv_table._get_numeric_data()
    e_nums = expr_table._get_numeric_data()
    # Pearson correlation coefficient (superfast)
    r = c_nums.corrwith(e_nums, axis=1)
    # Spearman, Kendall (slow)
    # NB: RuntimeWarning from numpy/scipy for divide-by-zero / inf / nan
    # ENH: with warnings.catch_warnings(): but just catch RuntimeWarning
    warnings.simplefilter('ignore', RuntimeWarning)
    rho, tau = zip(*[
        (spearmanr(cnv_row, expr_row, nan_policy='omit')[0],
         kendalltau(cnv_row, expr_row, nan_policy='omit')[0])
        for cnv_row, expr_row in zip(c_nums.values, e_nums.values)])

    result = pd.DataFrame({
        "kendall_t": tau,
        "pearson_r": r.to_numpy(),
        "spearman_r": rho,
    }, index=cnv_table.index).clip(lower=0).fillna(0)
    result.insert(0, 'hugo_gene',
            cnv_table['Hugo_Symbol'].fillna('').to_numpy())
    return result


def load_tcga_table(fname, shared_key):
    """Load TCGA expression/CNV data, keeping unique Entrez genes.

    Rows without an Entrez_Gene_Id value are dropped. Where a gene has multiple
    HUGO names but one Entrez_Gene_Id (i.e. multiple rows with the same
    Entrez_Gene_Id), only the sortest and then alphabetically first HUGO name is
    kept, ensuring Entrez_Gene_Id values are unique.
    """
    table = pd.read_csv(fname, sep='\t', dtype={shared_key: str}, na_filter=False)
    table = table[table[shared_key] != ''].astype({shared_key: int})
    before_pipe = before('|')
    sort_order = (table['Hugo_Symbol']
                  .apply(lambda x: (len(x), before_pipe(x)))
                  .argsort())
    table = (table.iloc[sort_order]
             .drop_duplicates(subset=shared_key)
             .set_index(shared_key)
             .sort_index(axis=0)
             .sort_index(axis=1))
    logging.info("Loaded %s with shape: %s", fname, table.shape)
    return table


def cnv_expression_correlate(args) -> None:
    table = correlate_cnv_expression(args.cnv_fname, args.expression_fname)
    table.to_csv(args.output or sys.stdout, sep='\t', index=True)
    if args.output:
        logging.info("Wrote %s with %s rows", args.output, len(table))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    arguments = argument_parsing()
    cnv_expression_correlate(arguments)


if __name__ == '__main__':
    main()
