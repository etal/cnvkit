"""Supporting functions for the 'reference' command."""
import collections
import logging
from concurrent import futures

import numpy as np
import pandas as pd
import pyfaidx
from skgenome import tabio, GenomicArray as GA

from . import core, fix, descriptives, params
from .cmdutil import read_cna
from .cnary import CopyNumArray as CNA


def do_reference_flat(targets, antitargets=None, fa_fname=None,
                      male_reference=False):
    """Compile a neutral-coverage reference from the given intervals.

    Combines the intervals, shifts chrX values if requested, and calculates GC
    and RepeatMasker content from the genome FASTA sequence.
    """
    ref_probes = bed2probes(targets)
    if antitargets:
        ref_probes.add(bed2probes(antitargets))
    # Set sex chromosomes by "reference" sex
    ref_probes['log2'] = ref_probes.expect_flat_log2(male_reference)
    ref_probes['depth'] = np.exp2(ref_probes['log2'])  # Shim
    # Calculate GC and RepeatMasker content for each probe's genomic region
    if fa_fname:
        gc, rmask = get_fasta_stats(ref_probes, fa_fname)
        ref_probes['gc'] = gc
        ref_probes['rmask'] = rmask
    else:
        logging.info("No FASTA reference genome provided; "
                     "skipping GC, RM calculations")
    ref_probes.sort_columns()
    return ref_probes


def bed2probes(bed_fname):
    """Create a neutral-coverage CopyNumArray from a file of regions."""
    regions = tabio.read_auto(bed_fname)
    table = regions.data.loc[:, ("chromosome", "start", "end")]
    table["gene"] = (regions.data["gene"] if "gene" in regions.data else '-')
    table["log2"] = 0.0
    table["spread"] = 0.0
    return CNA(table, {"sample_id": core.fbase(bed_fname)}) # todo: find all invocations of CNA() and forward param


def do_reference(target_fnames, antitarget_fnames=None, fa_fname=None,
                 male_reference=False, treat_par_on_chrx_as_autosomal_for_genome_build=None, female_samples=None,
                 do_gc=True, do_edge=True, do_rmask=True, do_cluster=False,
                 min_cluster_size=4):
    """Compile a coverage reference from the given files (normal samples)."""
    if antitarget_fnames:
        core.assert_equal("Unequal number of target and antitarget files given",
                          targets=len(target_fnames),
                          antitargets=len(antitarget_fnames))
    if not fa_fname:
        logging.info("No FASTA reference genome provided; "
                     "skipping GC, RM calculations")

    if female_samples is None:
        # NB: Antitargets are usually preferred for inferring sex, but might be
        # empty files, in which case no inference can be done. Since targets are
        # guaranteed to exist, infer from those first, then replace those
        # values where antitargets are suitable.
        sexes = infer_sexes(target_fnames, False, treat_par_on_chrx_as_autosomal_for_genome_build)
        if antitarget_fnames:
            a_sexes = infer_sexes(antitarget_fnames, False, treat_par_on_chrx_as_autosomal_for_genome_build)
            for sid, a_is_xx in a_sexes.items():
                t_is_xx = sexes.get(sid)
                if t_is_xx is None:
                    sexes[sid] = a_is_xx
                elif t_is_xx != a_is_xx and a_is_xx is not None:
                    logging.warning("Sample %s chromosomal X/Y ploidy looks "
                                    "like %s in targets but %s in antitargets; "
                                    "preferring antitargets",
                                    sid,
                                    "female" if t_is_xx else "male",
                                    "female" if a_is_xx else "male")
                    sexes[sid] = a_is_xx
    else:
        sexes = collections.defaultdict(lambda: female_samples)

    # TODO - refactor/inline this func here, once it works
    ref_probes = combine_probes(target_fnames, antitarget_fnames, fa_fname,
                                male_reference, treat_par_on_chrx_as_autosomal_for_genome_build,
                                sexes, do_gc, do_edge, do_rmask, do_cluster, min_cluster_size)
    warn_bad_bins(ref_probes)
    return ref_probes


def infer_sexes(cnn_fnames, is_haploid_x, treat_par_on_chrx_as_autosomal_for_genome_build):
    """Map sample IDs to inferred chromosomal sex, where possible.

    For samples where the source file is empty or does not include either sex
    chromosome, that sample ID will not be in the returned dictionary.
    """
    sexes = {}
    for fname in cnn_fnames:
        cnarr = read_cna(fname, treat_par_on_chrx_as_autosomal_for_genome_build)
        if cnarr:
            is_xx = cnarr.guess_xx(is_haploid_x)
            if is_xx is not None:
                sexes[cnarr.sample_id] = is_xx
    return sexes


def combine_probes(filenames, antitarget_fnames, fa_fname,
                   is_haploid_x, treat_par_on_chrx_as_autosomal_for_genome_build,
                   sexes, fix_gc, fix_edge, fix_rmask, do_cluster, min_cluster_size):
    """Calculate the median coverage of each bin across multiple samples.

    Parameters
    ----------
    filenames : list
        List of string filenames, corresponding to targetcoverage.cnn and
        antitargetcoverage.cnn files, as generated by 'coverage' or
        'import-picard'.
    fa_fname : str
        Reference genome sequence in FASTA format, used to extract GC and
        RepeatMasker content of each genomic bin.
    is_haploid_x : bool
    do_cluster : bool
    fix_gc : bool
    fix_edge : bool
    fix_rmask : bool

    Returns
    -------
    CopyNumArray
        One object summarizing the coverages of the input samples, including
        each bin's "average" coverage, "spread" of coverages, and GC content.
    """
    # Special parameters:
    # skip_low (when centering) = True for target; False for antitarget
    # do_edge  = as given for target; False for antitarget
    # do_rmask  = False for target; as given for antitarget
    ref_df, all_logr, all_depths = load_sample_block(
        filenames, fa_fname, is_haploid_x, treat_par_on_chrx_as_autosomal_for_genome_build,
        sexes, True, fix_gc, fix_edge, False)
    if antitarget_fnames:
        # XXX TODO ensure ordering matches targets!
        #   argsort on both -> same?
        anti_ref_df, anti_logr, anti_depths = load_sample_block(
            antitarget_fnames, fa_fname, is_haploid_x, treat_par_on_chrx_as_autosomal_for_genome_build,
            sexes, False, fix_gc, False, fix_rmask)
        ref_df = pd.concat([ref_df, anti_ref_df], ignore_index=True)
        all_logr = np.hstack([all_logr, anti_logr])
        all_depths = np.hstack([all_depths, anti_depths])

    stats_all = summarize_info(all_logr, all_depths)
    ref_df = ref_df.assign(**stats_all)

    if do_cluster:
        # Get extra cols, concat axis=1 here (DATAFRAME v-concat)
        sample_ids = [core.fbase(f) for f in filenames]
        if len(sample_ids) != len(all_logr) - 1:
            raise ValueError("Expected %d target coverage files (.cnn), got %d"
                             % (len(all_logr) - 1, len(sample_ids)))
        clustered_cols = create_clusters(all_logr, min_cluster_size, sample_ids)
        if clustered_cols:
            try:
                ref_df = ref_df.assign(**clustered_cols)
            except ValueError as exc:
                print("Reference:", len(ref_df.index))
                for cl_key, cl_col in clustered_cols.items():
                    print(cl_key, ":", len(cl_col))
                raise exc
        else:
            print("** Why weren't there any clustered cols?")

    ref_cna = CNA(ref_df, meta_dict={'sample_id': 'reference'})
    # NB: Up to this point, target vs. antitarget bins haven't been row-sorted.
    # Holding off until here ensures the cluster-specific log2 columns are
    # sorted with the same row order as the rest of the CNA dataframe.
    ref_cna.sort()
    ref_cna.sort_columns()
    # TODO figure out centering
    #ref_probes.center_all(skip_low=True)  # <-- on each cluster col, too?
    #   or just subtract the same amount from the cluster calls after figuring
    #   out the main one?
    return ref_cna


def load_sample_block(filenames, fa_fname,
                      is_haploid_x, treat_par_on_chrx_as_autosomal_for_genome_build,
                      sexes, skip_low, fix_gc, fix_edge, fix_rmask):
    r"""Load and summarize a pool of \*coverage.cnn files.

    Run separately for the on-target and (optional) antitarget bins.

    Returns
    -------
    ref_df : pandas.DataFrame
        All columns needed for the reference CNA object, including
        aggregate log2 and spread.
    all_logr : numpy.ndarray
        All sample log2 ratios, as a 2D matrix (rows=bins, columns=samples),
        to be used with do_cluster.
    """
    # Ensures samples' target and antitarget matrix columns are in the same
    # order, so they can be concatenated.
    # (We don't explicitly pair off each sample's target and antitarget .cnn
    # files; as long as the filename prefixes match, the resulting columns will
    # be consistent. Same is true for sample sex inference.)
    filenames = sorted(filenames, key=core.fbase)

    # Load coverage from target/antitarget files
    logging.info("Loading %s", filenames[0])
    cnarr1 = read_cna(filenames[0], treat_par_on_chrx_as_autosomal_for_genome_build)
    if not len(cnarr1):
        # Just create an empty array with the right columns
        col_names = ['chromosome', 'start', 'end', 'gene', 'log2', 'depth']
        if 'gc' in cnarr1 or fa_fname:
            col_names.append('gc')
        if fa_fname:
            col_names.append('rmask')
        col_names.append('spread')
        empty_df = pd.DataFrame.from_records([], columns=col_names)
        empty_logr = np.array([[]] * (len(filenames) + 1))
        empty_dp = np.array([[]] * len(filenames))
        return empty_df, empty_logr, empty_dp

    # Calculate GC and RepeatMasker content for each probe's genomic region
    ref_columns = {
        'chromosome': cnarr1.chromosome,
        'start': cnarr1.start,
        'end': cnarr1.end,
        'gene': cnarr1['gene'],
    }
    if fa_fname and (fix_rmask or fix_gc):
        gc, rmask = get_fasta_stats(cnarr1, fa_fname)
        if fix_gc:
            ref_columns['gc'] = gc
        if fix_rmask:
            ref_columns['rmask'] = rmask
    elif 'gc' in cnarr1 and fix_gc:
        # Reuse .cnn GC values if they're already stored (via import-picard)
        gc = cnarr1['gc']
        ref_columns['gc'] = gc

    # Make the sex-chromosome coverages of male and female samples compatible
    is_chr_x = cnarr1.chrx_filter()
    is_chr_y = (cnarr1.chromosome == cnarr1._chr_y_label)
    ref_flat_logr = cnarr1.expect_flat_log2(is_haploid_x)
    ref_edge_bias = fix.get_edge_bias(cnarr1, params.INSERT_SIZE)
    # Pseudocount of 1 "flat" sample
    all_depths = [cnarr1['depth'] if 'depth' in cnarr1
                  else np.exp2(cnarr1['log2'])]
    all_logr = [
        ref_flat_logr,
        bias_correct_logr(cnarr1, ref_columns, ref_edge_bias,
                          ref_flat_logr, sexes, is_chr_x, is_chr_y,
                          fix_gc, fix_edge, fix_rmask, skip_low)]

    # Load only coverage depths from the remaining samples

    procs = 1 # TODO: Add as param
    with futures.ProcessPoolExecutor(procs) as pool:
        args_iter = ((fname, treat_par_on_chrx_as_autosomal_for_genome_build, cnarr1, filenames[0], ref_columns, ref_edge_bias, ref_flat_logr, sexes, is_chr_x, is_chr_y, fix_gc, fix_edge, fix_rmask, skip_low) for fname in filenames[1:])
        for depths, logr in pool.map(_parallel_bias_correct_logr, args_iter):
            all_depths.append(depths)
            all_logr.append(logr)
    all_logr = np.vstack(all_logr)
    all_depths = np.vstack(all_depths)
    ref_df = pd.DataFrame.from_dict(ref_columns)
    return ref_df, all_logr, all_depths

def _parallel_bias_correct_logr(args):
    """Wrapper for parallel."""
    fname, treat_par_on_chrx_as_autosomal_for_genome_build, cnarr1, fname1, ref_columns, ref_edge_bias, ref_flat_logr, sexes, is_chr_x, is_chr_y, fix_gc, fix_edge, fix_rmask, skip_low = args
    logging.info("Loading %s", fname)
    cnarrx = read_cna(fname, treat_par_on_chrx_as_autosomal_for_genome_build)
    # Bin information should match across all files
    if not np.array_equal(
            cnarr1.data.loc[:, ('chromosome', 'start', 'end', 'gene')].values,
            cnarrx.data.loc[:, ('chromosome', 'start', 'end', 'gene')].values):
        raise RuntimeError("%s bins do not match those in %s"
                           % (fname, fname1))
    depths = cnarrx['depth'] if 'depth' in cnarrx else np.exp2(cnarrx['log2'])
    logr = bias_correct_logr(cnarrx, ref_columns, ref_edge_bias, ref_flat_logr,
                             sexes, is_chr_x, is_chr_y,
                             fix_gc, fix_edge, fix_rmask, skip_low)
    return (depths, logr)

def bias_correct_logr(cnarr, ref_columns, ref_edge_bias,
                      ref_flat_logr, sexes, is_chr_x, is_chr_y,
                      fix_gc, fix_edge, fix_rmask, skip_low):
    """Perform bias corrections on the sample."""
    cnarr.center_all(skip_low=skip_low)
    shift_sex_chroms(cnarr, sexes, ref_flat_logr, is_chr_x, is_chr_y)
    # Skip bias corrections if most bins have no coverage (e.g. user error)
    if (cnarr['log2'] > params.NULL_LOG2_COVERAGE - params.MIN_REF_COVERAGE
       ).sum() <= len(cnarr) // 2:
        logging.warning("WARNING: most bins have no or very low coverage; "
                        "check that the right BED file was used")
    else:
        # todo: output samlpe name while logging (due to parallelism)
        if 'gc' in ref_columns and fix_gc:
            logging.info("Correcting for GC bias...")
            cnarr = fix.center_by_window(cnarr, .1, ref_columns['gc'])
        if 'rmask' in ref_columns and fix_rmask:
            logging.info("Correcting for RepeatMasker bias...")
            cnarr = fix.center_by_window(cnarr, .1, ref_columns['rmask'])
        if fix_edge:
            logging.info("Correcting for density bias...")
            cnarr = fix.center_by_window(cnarr, .1, ref_edge_bias)
    return cnarr['log2']


def shift_sex_chroms(cnarr, sexes, ref_flat_logr, is_chr_x, is_chr_y):
    """Shift sample X and Y chromosomes to match the reference sex.

    Reference values::

        XY: chrX -1, chrY -1
        XX: chrX 0, chrY -1

    Plan::

        chrX:
        xx sample, xx ref: 0    (from 0)
        xx sample, xy ref: -= 1 (from -1)
        xy sample, xx ref: += 1 (from 0)    +1
        xy sample, xy ref: 0    (from -1)   +1
        chrY:
        xx sample, xx ref: = -1 (from -1)
        xx sample, xy ref: = -1 (from -1)
        xy sample, xx ref: 0    (from -1)   +1
        xy sample, xy ref: 0    (from -1)   +1

    """
    is_xx = sexes.get(cnarr.sample_id)
    cnarr['log2'] += ref_flat_logr
    if is_xx:
        # chrX has same ploidy as autosomes; chrY is just unusable noise
        cnarr[is_chr_y, 'log2'] = -1.0  # np.nan is worse
    else:
        # 1/2 #copies of each sex chromosome
        cnarr[is_chr_x | is_chr_y, 'log2'] += 1.0


def summarize_info(all_logr, all_depths):
    """Average & spread of log2ratios and depths for a group of samples.

    Can apply to all samples, or a given cluster of samples.
    """
    logging.info("Calculating average bin coverages")
    cvg_centers = np.apply_along_axis(descriptives.biweight_location, 0,
                                      all_logr)
    depth_centers = np.apply_along_axis(descriptives.biweight_location, 0,
                                        all_depths)
    logging.info("Calculating bin spreads")
    spreads = np.array([descriptives.biweight_midvariance(a, initial=i)
                        for a, i in zip(all_logr.T, cvg_centers)])

    result = {
        'log2': cvg_centers,
        'depth': depth_centers,
        'spread': spreads,
    }
    # TODO center the resulting log2
    #ref_df = pd.DataFrame.from_dict(ref_columns)
    #ref_cna = CNA.from_columns(ref_columns, {'sample_id': "reference"})
    return result


def create_clusters(logr_matrix, min_cluster_size, sample_ids):
    """Extract and summarize clusters of samples in logr_matrix.

    1. Calculate correlation coefficients between all samples (columns).
    2. Cluster the correlation matrix.
    3. For each resulting sample cluster (down to a minimum size threshold),
       calculate the central log2 value for each bin, similar to the full pool.
       Also print the sample IDs in each cluster, if feasible.

    Also recalculate and store the 'spread' of each cluster, though this might
    not be necessary/good.

    Return a DataFrame of just the log2 values. Column names are ``log2_i``
    where i=1,2,... .
    """
    from .cluster import markov, kmeans
    # Drop the pseudocount sample
    logr_matrix = logr_matrix[1:, :]
    print("Clustering", len(logr_matrix), "samples...")
    #clusters = markov(logr_matrix)
    clusters = kmeans(logr_matrix)
    cluster_cols = {}
    sample_ids = np.array(sample_ids)  # For easy indexing
    for i, clust_idx in enumerate(clusters):
        i += 1
        #print(len(clust_idx), clust_idx)
        if len(clust_idx) < min_cluster_size:
            logging.info("Skipping cluster #%d, size %d < min. %d",
                         i, len(clust_idx), min_cluster_size)
            continue
        logging.info("Summarizing cluster #%d of %d samples",
                        i, len(clust_idx))
        # List which samples are in each cluster
        samples = sample_ids[clust_idx]
        logging.info("\n".join(["\t" + s for s in samples]))
        # Calculate each cluster's summary stats
        clust_matrix = logr_matrix[clust_idx, :]
        # XXX re-add the pseudocount sample to each cluster? need benchmark
        clust_info = summarize_info(clust_matrix, [])
        cluster_cols.update({
            'log2_%d' % i: clust_info['log2'],
            'spread_%d' % i: clust_info['spread'],
        })
    return cluster_cols


def warn_bad_bins(cnarr, max_name_width=50):
    """Warn about target bins where coverage is poor.

    Prints a formatted table to stderr.
    """
    bad_bins = cnarr[fix.mask_bad_bins(cnarr)]
    fg_index = ~bad_bins['gene'].isin(params.ANTITARGET_ALIASES)
    fg_bad_bins = bad_bins[fg_index]
    if len(fg_bad_bins) > 0:
        bad_pct = (100 * len(fg_bad_bins)
                   / sum(~cnarr['gene'].isin(params.ANTITARGET_ALIASES)))
        logging.info("Targets: %d (%s) bins failed filters "
                     "(log2 < %s, log2 > %s, spread > %s)",
                     len(fg_bad_bins),
                     "%.4f" % bad_pct + '%',
                     params.MIN_REF_COVERAGE,
                     -params.MIN_REF_COVERAGE,
                     params.MAX_REF_SPREAD)
        if len(fg_bad_bins) < 500:
            gene_cols = min(max_name_width, max(map(len, fg_bad_bins['gene'])))
            labels = fg_bad_bins.labels()
            chrom_cols = max(labels.apply(len))
            last_gene = None
            for label, probe in zip(labels, fg_bad_bins):
                if probe.gene == last_gene:
                    gene = '  "'
                else:
                    gene = probe.gene
                    last_gene = gene
                if len(gene) > max_name_width:
                    gene = gene[:max_name_width-3] + '...'
                if 'rmask' in cnarr:
                    logging.info("  %s  %s  log2=%.3f  spread=%.3f  rmask=%.3f",
                                 gene.ljust(gene_cols),
                                 label.ljust(chrom_cols),
                                 probe.log2, probe.spread, probe.rmask)
                else:
                    logging.info("  %s  %s  log2=%.3f  spread=%.3f",
                                 gene.ljust(gene_cols),
                                 label.ljust(chrom_cols),
                                 probe.log2, probe.spread)

    # Count the number of BG bins dropped, too (names are all "Antitarget")
    bg_bad_bins = bad_bins[~fg_index]
    if len(bg_bad_bins) > 0:
        bad_pct = (100 * len(bg_bad_bins)
                   / sum(cnarr['gene'].isin(params.ANTITARGET_ALIASES)))
        logging.info("Antitargets: %d (%s) bins failed filters",
                     len(bg_bad_bins), "%.4f" % bad_pct + '%')


def get_fasta_stats(cnarr, fa_fname):
    """Calculate GC and RepeatMasker content of each bin in the FASTA genome."""
    logging.info("Calculating GC and RepeatMasker content in %s ...", fa_fname)
    gc_rm_vals = [calculate_gc_lo(subseq)
                  for subseq in fasta_extract_regions(fa_fname, cnarr)]
    gc_vals, rm_vals = zip(*gc_rm_vals)
    return np.asfarray(gc_vals), np.asfarray(rm_vals)


def calculate_gc_lo(subseq):
    """Calculate the GC and lowercase (RepeatMasked) content of a string."""
    cnt_at_lo = subseq.count('a') + subseq.count('t')
    cnt_at_up = subseq.count('A') + subseq.count('T')
    cnt_gc_lo = subseq.count('g') + subseq.count('c')
    cnt_gc_up = subseq.count('G') + subseq.count('C')
    tot = float(cnt_gc_up + cnt_gc_lo + cnt_at_up + cnt_at_lo)
    if not tot:
        return 0.0, 0.0
    frac_gc = (cnt_gc_lo + cnt_gc_up) / tot
    frac_lo = (cnt_at_lo + cnt_gc_lo) / tot
    return frac_gc, frac_lo


def fasta_extract_regions(fa_fname, intervals):
    """Extract an iterable of regions from an indexed FASTA file.

    Input: FASTA file name; iterable of (seq_id, start, end) (1-based)
    Output: iterable of string sequences.
    """
    with pyfaidx.Fasta(fa_fname, as_raw=True) as fa_file:
        for chrom, subarr in intervals.by_chromosome():
            logging.info("Extracting sequences from chromosome %s", chrom)
            for _chrom, start, end in subarr.coords():
                yield fa_file[_chrom][int(start):int(end)]


def reference2regions(refarr):
    """Split reference into target and antitarget regions."""
    is_bg = (refarr['gene'].isin(params.ANTITARGET_ALIASES))
    regions = GA(refarr.data.loc[:, ('chromosome', 'start', 'end', 'gene')],
                 {'sample_id': 'reference'})
    targets = regions[~is_bg]
    antitargets = regions[is_bg]
    return targets, antitargets
