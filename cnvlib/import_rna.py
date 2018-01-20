from __future__ import absolute_import, division, print_function
import logging
import os

import numpy as np
import pandas as pd

from . import rna


def do_import_rna(gene_count_fnames, in_format, gene_resource_fname,
                  correlations_fname=None, normal_fnames=()):
    """Convert a cohort of per-gene read counts to CNVkit .cnr format.

    The expected data source is TCGA gene-level expression counts for individual
    samples, but other sources should be fine, too.
    """
    # Deduplicate and ensure all normals are included in the analysis
    gene_count_fnames = sorted(set(gene_count_fnames) + set(normal_fnames))

    if in_format == 'rsem':
        sample_counts, tx_lengths = aggregate_rsem(gene_count_fnames)
    elif in_format == 'counts':
        sample_counts = aggregate_gene_counts(gene_count_fnames)
        tx_lengths = None
    else:
        raise RuntimeError("Unrecognized input format name: %r" % in_format)
    sample_counts = rna.filter_probes(sample_counts)

    logging.info("Loading gene metadata" +
                 (" and TCGA gene expression/CNV profiles"
                  if correlations_fname else ""))
    gene_info = rna.load_gene_info(gene_resource_fname, correlations_fname)

    logging.info("Aligning gene info to sample gene counts")
    normal_ids = [os.path.basename(f).split('.')[0] for f in normal_fnames]
    gene_info, sample_counts, sample_data_log2 = rna.align_gene_info_to_samples(
        gene_info, sample_counts, tx_lengths, normal_ids)

    # Summary table has log2-normalized values, not raw counts
    # ENH show both, with column header suffixes to distinguish?
    all_data = pd.concat([gene_info, sample_data_log2], axis=1)
    # CNVkit files have both absolute and log2-normalized read counts
    cnrs = rna.attach_gene_info_to_cnr(sample_counts, sample_data_log2,
                                       gene_info)
    cnrs = (rna.correct_cnr(cnr) for cnr in cnrs)
    return all_data, cnrs


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


def aggregate_rsem(fnames):
    """Pull out the expected read counts from each RSEM file.

    The format of RSEM's ``*_rsem.genes.results`` output files is tab-delimited
    with a header row. We extract the Ensembl gene ID, expected read counts, and
    transcript lengths from each file.

    Returns
    -------
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
