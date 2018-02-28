"""RNA expression levels quantified per gene.

Process per-gene expression levels, or the equivalent, by cohort.
"""
from __future__ import absolute_import, division, print_function
from builtins import zip

import logging
import sys

import numpy as np
import pandas as pd
from scipy.stats import gmean

from .cnary import CopyNumArray as CNA
from .fix import center_by_window
# from .descriptives import biweight_midvariance
# from .params import NULL_LOG2_COVERAGE

NULL_LOG2_COVERAGE = -5


def before(char):
    def selector(string):
        return string.split(char, 1)[0]
    return selector


def filter_probes(sample_counts):
    """Filter probes to only include high-quality, transcribed genes.

    The human genome has ~25,000 protein coding genes, yet the RSEM output
    includes ~58,000 rows. Some of these rows correspond to read counts over
    control probes (e.g.  spike-in sequences). Some rows correspond to poorly
    mapped genes in contigs that have not been linked to the 24 chromosomes
    (e.g. HLA region). Others correspond to pseudo-genes and non-coding genes.
    For the purposes of copy number inference, these rows are best removed.
    """
    gene_medians = sample_counts.median(axis=1)
    ## Make sure the gene has detectable transcript in at least half of samples
    is_mostly_transcribed = (gene_medians >= 1.0)
    print("Dropping", (~is_mostly_transcribed).sum(), "/",
          len(is_mostly_transcribed), "rarely expressed genes from input samples")
    return sample_counts[is_mostly_transcribed]


def load_gene_info(gene_resource, corr_fname, default_r=.1):
    """Read gene info from BioMart, and optionally TCGA, into a dataframe.

    RSEM only outputs the Ensembl ID. We have used the BioMART tool in Ensembl
    to export a list of all Ensembl genes with a RefSeq mRNA (meaning it is high
    quality, curated, and bona fide gene) and resides on chromosomes 1-22, X, or
    Y. The tool also outputs the GC content of the gene, chromosomal coordinates
    of the gene, and HUGO gene symbol.

    The gene resource input can be obtained from a resource bundle we provide
    (for reference genome hg19) or generated from BioMart.
    """
    # Load the gene info file and clean up column names
    # Original columns:
    # "Gene stable ID" -- Ensembl ID
    # "Gene % GC content"
    # "Chromosome/scaffold name"
    # "Gene start (bp)"
    # "Gene end (bp)"
    # "Gene name"
    # "NCBI gene ID"
    # "Transcript length (including UTRs and CDS)"
    # "Transcript support level (TSL)"
    info_col_names = ['gene_id', 'gc', 'chromosome', 'start', 'end',
                      'gene', 'entrez_id', 'tx_length', 'tx_support']
    gene_info = (pd.read_table(gene_resource, names=info_col_names, header=1,
                               converters={'gene_id': before('.'),
                                           'tx_support': tsl2int,
                                           'gc': lambda x: float(x)/100})
                 .sort_values('gene_id'))
    print("Loaded", gene_resource, "with shape:", gene_info.shape,
          file=sys.stderr)

    if corr_fname:
        corr_table = load_cnv_expression_corr(corr_fname)

        gi_corr = gene_info.join(corr_table, on='entrez_id', how='left')
        if not gi_corr['gene_id'].is_unique:
            unique_idx = gi_corr.groupby('gene_id').apply(dedupe_ens_hugo)
            gi_corr = gi_corr.loc[unique_idx]

        if not gene_info['entrez_id'].is_unique:
            # Treat untraceable entrez_id duplicates (scarce at this point, ~15)
            # as unknown correlation, i.e. default, as if not in TCGA
            gi_corr.loc[tuple(locate_entrez_dupes(gi_corr)),
                        ('pearson_r', 'spearman_r', 'kendall_t')] = default_r

        # Genes w/o TCGA info get default correlation 0.1 (= .5*corr.median())
        gene_info = gi_corr.fillna({'pearson_r': default_r,
                                    'spearman_r': default_r,
                                    'kendall_t': default_r,
                                   })

    elif not gene_info['gene_id'].is_unique:
        # Simpler deduplication, not using TCGA fields
        unique_idx = gene_info.groupby('gene_id').apply(dedupe_ens_no_hugo)
        gene_info = gene_info.loc[unique_idx]

    print("Trimmed gene info table to shape:", gene_info.shape)
    assert gene_info['gene_id'].is_unique
    gene_info['entrez_id'] = gene_info['entrez_id'].fillna(0).astype('int')
    return gene_info.set_index('gene_id')


def load_cnv_expression_corr(fname):
    shared_key = 'Entrez_Gene_Id'
    table = (pd.read_table(fname, dtype={shared_key: str}, na_filter=False)
             .set_index(shared_key))
    print("Loaded", fname, "with shape:", table.shape, file=sys.stderr)
    return table


def tsl2int(tsl):
    """Convert an Ensembl Transcript Support Level (TSL) code to an integer.

    The code has the format "tsl([1-5]|NA)".

    See: https://www.ensembl.org/Help/Glossary?id=492
    """
    if tsl in (np.nan, "", "tslNA"):
        return 0
    # Remove "tsl" prefix
    assert tsl[:3] == "tsl"
    value = tsl[3:]
    # Remove possible suffix " (assigned to previous version _)"
    if len(value) > 2:
        value = value[:2].rstrip()
    if value == "NA":
        return 0
    return int(value)


def dedupe_ens_hugo(dframe):
    """Emit the "best" index from a group of the same Ensembl ID.

    The RSEM gene rows are the data of interest, and they're associated with
    Ensembl IDs to indicate the transcribed gene being measured in each row.
    The BioMart table of gene info can have duplicate rows for Ensembl ID, which
    would cause duplicate rows in RSEM import if joined naively.  So, we need to
    select a single row for each group of "gene info" rows with the same Ensembl
    ID (column 'gene_id').

    The keys we can use for this are:

    - Entrez ID ('entrez_id')
    - Ensembl gene name ('gene')
    - Entrez gene name ('hugo_gene')

    Entrez vs. Ensembl IDs and gene names are potentially many-to-many, e.g.
    CALM1/2/3. However, if we also require that the Ensembl and HUGO gene names
    match within a group, that (usually? always?) results in a unique row
    remaining.

    (Example: CALM1/2/3 IDs are many-to-many, but of the 3 Entrez IDs associated
    with Ensembl's CALM1, only 1 is called CALM1 in the Entrez/corr. table.)

    Failing that (no matches or multiple matches), prefer a lower Entrez ID,
    because well-characterized, protein-coding genes tend to have been
    discovered and accessioned first.
    """
    if len(dframe) == 1:
        return dframe.index[0]
    match_gene = dframe[dframe['gene'] == dframe['hugo_gene']]
    if len(match_gene) == 1:
        return match_gene.index[0]
    if len(match_gene) > 1:
        # Take lowest Entrez ID of the matched
        return dedupe_tx(match_gene)
    # No matching names -> lowest Entrez ID key
    return dedupe_tx(dframe)


def dedupe_ens_no_hugo(dframe):
    """Deduplicate Ensembl ID using Entrez ID but not HUGO gene name."""
    if len(dframe) == 1:
        return dframe.index[0]
    return dedupe_tx(dframe)


def dedupe_tx(dframe):
    """Deduplicate table rows to select one transcript length per gene.

    Choose the lowest-number Entrez ID and the transcript with the greatest
    support (primarily) and length (secondarily).

    This is done at the end of Ensembl ID deduplication, after filtering on gene
    names and for single-row tables.

    Returns an integer row index corresponding to the original table.
    """
    # NB: Transcripts are many-to-1 vs. Ensembl ID, not Entrez ID, so each
    # unique Entrez ID here will have the same collection of transcripts
    # associated with it.
    return (dframe.sort_values(['entrez_id', 'tx_support', 'tx_length'],
                               ascending=[True, False, False],
                               na_position='last')
            .index[0])


def locate_entrez_dupes(dframe):
    """In case the same Entrez ID was assigned to multiple Ensembl IDs.

    Use HUGO vs. HGNC name again, similar to `dedupe_hugo`, but instead of
    emiting the indices of the rows to keep, emit the indices of the extra rows
    -- their correlation values will then be filled in with a default value
    (np.nan or 0.1).

    It will then be as if those genes hadn't appeared in the TCGA tables at all,
    i.e. CNV-expression correlation is unknown, but all entries are still
    retained in the BioMart table (gene_info).
    """
    for _key, group in dframe.groupby('entrez_id'):
        if len(group) == 1:
            continue
        match_gene_idx = (group['gene'] == group['hugo_gene'])
        match_gene_cnt = match_gene_idx.sum()
        if match_gene_cnt == 1:
            for mismatch_idx in group.index[~match_gene_idx]:
                yield mismatch_idx
        else:
            # Keep the lowest Ensemble ID (of the matched, if any)
            if match_gene_cnt: # >1
                keepable = group[match_gene_idx]
            else:
                keepable = group
            idx_to_keep = keepable.sort_values('gene_id').index.values[0]
            for idx in group.index:
                if idx != idx_to_keep:
                    yield idx


def align_gene_info_to_samples(gene_info, sample_counts, tx_lengths,
                               normal_ids):
    """Align columns and sort.

    Also calculate weights and add to gene_info as 'weight', along with
    transcript lengths as 'tx_length'.
    """
    print("Dimensions: gene_info=%s, sample_counts=%s"
          % (gene_info.shape, sample_counts.shape))
    sc, gi = sample_counts.align(gene_info, join='inner', axis=0)
    gi = gi.sort_values(by=['chromosome', 'start'])
    sc = sc.loc[gi.index]

    if tx_lengths is not None:
        # Replace the existing tx_lengths from gene_resource
        # (RSEM results have this, TCGA gene counts don't)
        gi['tx_length'] = tx_lengths.loc[gi.index]

    # Calculate per-gene weights similarly to cnvlib.fix
    # NB: chrX doesn't need special handling because with X-inactivation,
    # expression should be similar in male and female samples, i.e. neutral is 0
    logging.info("Weighting genes with below-average read counts")
    gene_counts = sc.median(axis=1)
    weights = [np.sqrt((gene_counts / gene_counts.quantile(.75)).clip_upper(1))]

    logging.info("Calculating normalized gene read depths")
    sample_depths_log2 = normalize_read_depths(sc.divide(gi['tx_length'],
                                                         axis=0),
                                               normal_ids)

    logging.info("Weighting genes by spread of read depths")
    gene_spreads = sample_depths_log2.std(axis=1)
    weights.append(gene_spreads)

    corr_weights = []
    for corr_col in ('spearman_r', 'pearson_r', 'kendall_t'):
        if corr_col in gi:
            logging.info("Weighting genes by %s correlation coefficient",
                         corr_col)
            corr_weights.append(gi[corr_col].values)
    if corr_weights:
        weights.append(np.vstack(corr_weights).mean(axis=0))

    weight = gmean(np.vstack(weights), axis=0)
    gi['weight'] = weight / weight.max()
    if gi['weight'].isnull().all():
        gi['weight'] = 1.0
    logging.debug(" --> final zeros: %d / %d",
                  (gi['weight'] == 0).sum(), len(gi))
    return gi, sc, sample_depths_log2


def normalize_read_depths(sample_depths, normal_ids):
    """Normalize read depths within each sample.

    Some samples have higher sequencing depth and therefore read depths need to
    be normalized within each sample. TCGA recommends an upper quartile
    normalization.

    After normalizing read depths within each sample, normalize (median-center)
    within each gene, across samples.

    Finally, convert to log2 ratios.
    """
    # TODO use normal_ids as a control here
    #   e.g. subtract their IQR after the loop?
    assert sample_depths.values.sum() > 0
    sample_depths = sample_depths.fillna(0)
    for _i in range(4):
        # By-sample: 75%ile  among all genes
        q3 = sample_depths.quantile(.75)
        sample_depths /= q3
        # By-gene: median among all samples
        sm = sample_depths.median(axis=1)
        sample_depths = sample_depths.divide(sm, axis=0)
        # print("  round", _i, "grand mean:", sample_depths.values.mean(),
        #       ": sample-wise denom =", q3.mean(),
        #       ", gene-wise denom =", sm.mean())
    # Finally, convert normalized read depths to log2 scale
    return safe_log2(sample_depths, NULL_LOG2_COVERAGE)


def safe_log2(values, min_log2):
    """Transform values to log2 scale, safely handling zeros.

    Parameters
    ----------
    values : np.array
        Absolute-scale values to transform. Should be non-negative.
    min_log2 : float
        Assign input zeros this log2-scaled value instead of -inf. Rather than
        hard-clipping, input values near 0 (especially below 2^min_log2) will be
        squeezed a bit above `min_log2` in the log2-scale output.
    """
    absolute_shift = 2 ** min_log2
    return np.log2(values + absolute_shift)


def attach_gene_info_to_cnr(sample_counts, sample_data_log2, gene_info,
                            read_len=100):
    """Join gene info to each sample's log2 expression ratios.

    Add the Ensembl gene info to the aggregated gene expected read counts,
    dropping genes that are not in the Ensembl table
    I.e., filter probes down to those genes that have names/IDs in the gene
    resource table.

    Split out samples to individual .cnr files, keeping (most) gene info.
    """
    gi_cols = ['chromosome', 'start', 'end', 'gene', 'gc', 'tx_length', 'weight']
    cnr_info = gene_info.loc[:, gi_cols]
    # Fill NA fields with the lowest finite value in the same row.
    # Only use NULL_LOG2_COVERAGE if all samples are NA / zero-depth.
    gene_minima = sample_data_log2.min(axis=1, skipna=True)
    assert not gene_minima.hasnans
    for (sample_id, sample_col), (_sid_log2, sample_log2) \
            in zip(sample_counts.iteritems(), sample_data_log2.iteritems()):
        tx_len = cnr_info.tx_length
        sample_depth = (read_len * sample_col / tx_len).rename("depth")
        sample_log2 = sample_log2.fillna(gene_minima).rename("log2")
        cdata = (pd.concat([cnr_info, sample_depth, sample_log2], axis=1)
                 .reset_index(drop=True))
        cnr = CNA(cdata, {'sample_id': sample_id})
        cnr.sort_columns()
        yield cnr


def correct_cnr(cnr):
    """Apply bias corrections & smoothing.

    - Biases: 'gc', 'length'
    - Smoothing: rolling triangle window using weights.
    """
    cnr.center_all()
    # Biases, similar to stock CNVkit
    if 'gc' in cnr:
        cnr = center_by_window(cnr, .1, cnr['gc'])
    if 'tx_length' in cnr:
        cnr = center_by_window(cnr, .1, cnr['tx_length'])
    cnr.center_all()
    return cnr
