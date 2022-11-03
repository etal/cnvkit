"""Supporting functions for the 'fix' command."""
import logging

import numpy as np
import pandas as pd

from . import descriptives, params, smoothing


def do_fix(target_raw, antitarget_raw, reference,
           do_gc=True, do_edge=True, do_rmask=True, do_cluster=False):
    """Combine target and antitarget coverages and correct for biases."""
    # Load, recenter and GC-correct target & antitarget probes separately
    logging.info("Processing target: %s", target_raw.sample_id)
    cnarr, ref_matched = load_adjust_coverages(target_raw, reference,
                                               True, do_gc, do_edge, False)
    logging.info("Processing antitarget: %s", antitarget_raw.sample_id)
    anti_cnarr, ref_anti = load_adjust_coverages(antitarget_raw, reference,
                                                 False, do_gc, False, do_rmask)
    if len(anti_cnarr):
        # Combine target and antitarget bins
        cnarr.add(anti_cnarr)
        ref_matched.add(ref_anti)

    # Find reference clusters, if requested
    log2_key = 'log2'
    spread_key = 'spread'
    if do_cluster:
        ref_log2_cols = [col for col in ref_matched.data.columns
                         if col == "log2"
                         or col.startswith("log2")]
        if len(ref_log2_cols) == 1:
            logging.info("Reference does not contain any sub-clusters; "
                         "using %s", log2_key)
        else:
            # Get correlations between test sample and each reference cluster
            corr_coefs = np.array([cnarr.log2.corr(ref_matched[ref_col])
                                   for ref_col in ref_log2_cols])
            ordered = [(k, r) for r, k in sorted(zip(corr_coefs, ref_log2_cols),
                                                 reverse=True)]
            logging.info("Correlations with each cluster:\n\t%s",
                         "\n\t".join(["{}\t: {}".format(k, r)
                                      for k, r in ordered]))
            log2_key = ordered[0][0]
            if log2_key.startswith('log2_'):
                suffix = log2_key.split('_', 1)[1]
                spread_key = 'spread_' + suffix
            logging.info(" -> Choosing columns %r and %r", log2_key, spread_key)

    # Normalize coverages according to the reference
    # (Subtract the reference log2 copy number to get the log2 ratio)
    cnarr.data['log2'] -= ref_matched[log2_key]
    cnarr = apply_weights(cnarr, ref_matched, log2_key, spread_key)
    cnarr.center_all(skip_low=True)
    return cnarr


def load_adjust_coverages(cnarr, ref_cnarr, skip_low,
                          fix_gc, fix_edge, fix_rmask):
    """Load and filter probe coverages; correct using reference and GC."""
    if 'gc' in cnarr:
        # Don't choke on Picard-derived files that have the GC column
        cnarr = cnarr.keep_columns(cnarr._required_columns + ('depth',))

    # No corrections needed if there are no data rows (e.g. no antitargets)
    if not len(cnarr):
        return cnarr, ref_cnarr[:0]

    ref_matched = match_ref_to_sample(ref_cnarr, cnarr)

    # Drop bins that had poor coverage in the pooled reference
    ok_cvg_indices = ~mask_bad_bins(ref_matched)
    logging.info("Keeping %d of %d bins", sum(ok_cvg_indices), len(ref_matched))
    cnarr = cnarr[ok_cvg_indices]
    ref_matched = ref_matched[ok_cvg_indices]

    # Apply corrections for known systematic biases in coverage
    cnarr.center_all(skip_low=skip_low)
    # Skip bias corrections if most bins have no coverage (e.g. user error)
    if (cnarr['log2'] > params.NULL_LOG2_COVERAGE - params.MIN_REF_COVERAGE
        ).sum() <= len(cnarr) // 2:
        logging.warning("WARNING: most bins have no or very low coverage; "
                        "check that the right BED file was used")
    else:
        cnarr_index_reset = False
        if fix_gc:
            if 'gc' in ref_matched:
                logging.info("Correcting for GC bias...")
                cnarr = center_by_window(cnarr, .1, ref_matched['gc'])
                cnarr_index_reset = True
            else:
                logging.warning("WARNING: Skipping correction for GC bias")
        if fix_edge:
            logging.info("Correcting for density bias...")
            edge_bias = get_edge_bias(cnarr, params.INSERT_SIZE)
            cnarr = center_by_window(cnarr, .1, edge_bias)
            cnarr_index_reset = True
        if fix_rmask:
            if 'rmask' in ref_matched:
                logging.info("Correcting for RepeatMasker bias...")
                cnarr = center_by_window(cnarr, .1, ref_matched['rmask'])
                cnarr_index_reset = True
            else:
                logging.warning("WARNING: Skipping correction for "
                                "RepeatMasker bias")
        if cnarr_index_reset:
            ref_matched.data.reset_index(drop=True, inplace=True)
    return cnarr, ref_matched


def mask_bad_bins(cnarr):
    """Flag the bins with excessively low or inconsistent coverage.

    Returns
    -------
    np.array
        A boolean array where True indicates bins that failed the checks.
    """
    mask = ((cnarr['log2'] < params.MIN_REF_COVERAGE) |
            (cnarr['log2'] > -params.MIN_REF_COVERAGE) |
            (cnarr['spread'] > params.MAX_REF_SPREAD))
    if 'depth' in cnarr:
        mask |= cnarr['depth'] == 0
    if 'gc' in cnarr:
        assert (params.GC_MIN_FRACTION >= 0 and params.GC_MIN_FRACTION <= 1)
        assert (params.GC_MAX_FRACTION >= 0 and params.GC_MAX_FRACTION <= 1)
        lower_gc_bound = min(params.GC_MIN_FRACTION, params.GC_MAX_FRACTION)
        upper_gc_bound = max(params.GC_MIN_FRACTION, params.GC_MAX_FRACTION)
        mask |= (cnarr['gc'] > upper_gc_bound) | (cnarr['gc'] < lower_gc_bound)
    return mask


def match_ref_to_sample(ref_cnarr, samp_cnarr):
    """Filter the reference bins to match the sample (target or antitarget)."""
    # Assign each bin a unique string ID based on genomic coordinates
    samp_labeled = samp_cnarr.data.set_index(pd.Index(samp_cnarr.coords()))
    ref_labeled = ref_cnarr.data.set_index(pd.Index(ref_cnarr.coords()))
    for dset, name in ((samp_labeled, "sample"),
                       (ref_labeled, "reference")):
        dupes = dset.index.duplicated()
        if dupes.any():
            raise ValueError(("Duplicated genomic coordinates in {} set. Total duplicated regions: {}, starting with:\n"
                              "{}.").format(name, len(dset.index[dupes]), "\n".join(map(str, dset.index[dupes][:10]))))
    # Take the reference bins with IDs identical to those in the sample
    ref_matched = ref_labeled.reindex(index=samp_labeled.index)
    # Check for signs that the wrong reference was used
    num_missing = pd.isnull(ref_matched.start).sum()
    if num_missing > 0:
        raise ValueError("Reference is missing %d bins found in %s"
                         % (num_missing, samp_cnarr.sample_id))
    x = ref_cnarr.as_dataframe(ref_matched.reset_index(drop=True)
                               .set_index(samp_cnarr.data.index))
    return x


def center_by_window(cnarr, fraction, sort_key):
    """Smooth out biases according to the trait specified by sort_key.

    E.g. correct GC-biased bins by windowed averaging across similar-GC
    bins; or for similar interval sizes.
    """
    # Separate neighboring bins that could have the same key
    # (to avoid re-centering actual CNV regions -- only want an independently
    # sampled subset of presumably overall-CN-neutral bins)
    df = cnarr.data.reset_index(drop=True)
    np.random.seed(0xA5EED)
    shuffle_order = np.random.permutation(df.index)
    #df = df.reindex(shuffle_order)
    df = df.iloc[shuffle_order]
    # Apply the same shuffling to the key array as to the target probe set
    if isinstance(sort_key, pd.Series):
        # XXX why
        sort_key = sort_key.values
    sort_key = sort_key[shuffle_order]
    # Sort the data according to the specified parameter
    order = np.argsort(sort_key, kind='mergesort')
    df = df.iloc[order]
    biases = smoothing.rolling_median(df['log2'], fraction)
    # biases = smoothing.savgol(df['log2'], fraction)
    df['log2'] -= biases
    fixarr = cnarr.as_dataframe(df)
    fixarr.sort()
    return fixarr


def get_edge_bias(cnarr, margin):
    """Quantify the "edge effect" of the target tile and its neighbors.

    The result is proportional to the change in the target's coverage due to
    these edge effects, i.e. the expected loss of coverage near the target
    edges and, if there are close neighboring tiles, gain of coverage due
    to "spill over" reads from the neighbor tiles.

    (This is not the actual change in coverage. This is just a tribute.)
    """
    output_by_chrom = []
    for _chrom, subarr in cnarr.by_chromosome():
        tile_starts = subarr['start'].values
        tile_ends = subarr['end'].values
        tgt_sizes = tile_ends - tile_starts
        # Calculate coverage loss at (both edges of) each tile
        losses = edge_losses(tgt_sizes, margin)
        # Find tiled intervals within a margin (+/- bp) of the given probe
        # (excluding the probe itself), then calculate the relative coverage
        # "gain" due to the neighbors, if any
        gap_sizes = tile_starts[1:] - tile_ends[:-1]
        ok_gaps_mask = (gap_sizes < margin)
        ok_gaps = gap_sizes[ok_gaps_mask]
        left_gains = edge_gains(tgt_sizes[1:][ok_gaps_mask], ok_gaps, margin)
        right_gains = edge_gains(tgt_sizes[:-1][ok_gaps_mask], ok_gaps, margin)
        gains = np.zeros(len(subarr))
        gains[np.concatenate([[False], ok_gaps_mask])] += left_gains
        gains[np.concatenate([ok_gaps_mask, [False]])] += right_gains
        output_by_chrom.append(gains - losses)
    return pd.Series(np.concatenate(output_by_chrom), index=cnarr.data.index)


def edge_losses(target_sizes, insert_size):
    """Calculate coverage losses at the edges of baited regions.

    Letting i = insert size and t = target size, the proportional loss of
    coverage near the two edges of the baited region (combined) is:

    .. math :: i/2t

    If the "shoulders" extend outside the bait $(t < i), reduce by:

    .. math :: (i-t)^2 / 4it

    on each side, or (i-t)^2 / 2it total.
    """
    losses = insert_size / (2 * target_sizes)
    # Drop the shoulder part that would extend past the bait
    small_mask = (target_sizes < insert_size)
    t_small = target_sizes[small_mask]
    losses[small_mask] -= ((insert_size - t_small)**2
                           / (2 * insert_size * t_small))
    return losses


def edge_gains(target_sizes, gap_sizes, insert_size):
    """Calculate coverage gain from neighboring baits' flanking reads.

    Letting i = insert size, t = target size, g = gap to neighboring bait,
    the gain of coverage due to a nearby bait, if g < i, is::

    .. math :: (i-g)^2 / 4it

    If the neighbor flank extends beyond the target (t+g < i), reduce by::

    .. math :: (i-t-g)^2 / 4it

    If a neighbor overlaps the target, treat it as adjacent (gap size 0).
    """
    if not (gap_sizes <= insert_size).all():
        raise ValueError("Gaps greater than insert size:\n" +
                         gap_sizes[gap_sizes > insert_size].head())
    gap_sizes = np.maximum(0, gap_sizes)
    gains = ((insert_size - gap_sizes)**2
             / (4 * insert_size * target_sizes))
    # Drop the flank part that extends past this baited region
    past_other_side_mask = (target_sizes + gap_sizes < insert_size)
    g_past = gap_sizes[past_other_side_mask]
    t_past = target_sizes[past_other_side_mask]
    gains[past_other_side_mask] -= ((insert_size - t_past - g_past)**2
                                    / (4 * insert_size * t_past))
    return gains


def apply_weights(cnarr, ref_matched, log2_key, spread_key, epsilon=1e-4):
    """Calculate weights for each bin.

    Bin weight is an estimate of (1 - variance) and within the range
    ``(0, 1]``.

    Weights are derived from:

    - Each bin's size
    - Sample's genome-wide average (on/off-target) coverage depth
    - Sample's genome-wide observed (on/off-target) bin variances

    And with a pooled reference:

    - Each bin's coverage depth in the reference
    - The "spread" column of the reference (approx. stdev)

    These estimates of variance assume the number of aligned reads per bin
    follows a Poisson distribution, approximately log-normal.

    Parameters
    ----------
    cnarr : CopyNumArray
        Sample bins.
    ref_match : CopyNumArray
        Reference bins.
    log2_key : string
        The 'log2' column name in the reference to use. A clustered reference
        may have a suffix indicating the cluster, e.g. "log2_1".
    spread_key : string
        The 'spread' or 'spread_<cluster_id>' column name to use.
    epsilon : float
        Minimum value for bin weights, to avoid 0-weight bins causing errors
        later during segmentation. (CBS doesn't allow 0-weight bins.)

    Returns: The input `cnarr` with a `weight` column added.
    """
    # Weight by sample-level features -- works for flat reference, too
    logging.debug("Weighting bins by size and overall variance in sample")
    simple_wt = np.zeros(len(cnarr))
    # Calculate separately for on-, off-target bins
    is_anti = cnarr['gene'].isin(params.ANTITARGET_ALIASES)
    tgt_cna = cnarr[~is_anti]
    tgt_var = descriptives.biweight_midvariance(tgt_cna
                                                .drop_low_coverage()
                                                .residuals()) ** 2
    bin_sz = np.sqrt(tgt_cna['end'] - tgt_cna['start'])
    tgt_simple_wts = 1 - tgt_var / (bin_sz / bin_sz.mean())
    simple_wt[~is_anti] = tgt_simple_wts

    if is_anti.any():
        # Check for a common user error
        anti_cna = cnarr[is_anti]
        anti_ok = anti_cna.drop_low_coverage()
        frac_anti_low = 1 - (len(anti_ok) / len(anti_cna))
        if frac_anti_low > .5:
            # Off-target bins are mostly garbage -- skip reweighting
            logging.warning("WARNING: Most antitarget bins ({:.2f}%, {:d}/{:d})"
                            " have low or no coverage; is this amplicon/WGS?"
                            .format(100 * frac_anti_low,
                                    len(anti_cna) - len(anti_ok),
                                    len(anti_cna)))

        anti_var = descriptives.biweight_midvariance(anti_ok.residuals()) ** 2
        anti_bin_sz = np.sqrt(anti_cna['end'] - anti_cna['start'])
        anti_simple_wts = 1 - anti_var / (anti_bin_sz / anti_bin_sz.mean())
        simple_wt[is_anti] = anti_simple_wts

        # Report any difference in bin set variability
        var_ratio = max(tgt_var, .01) / max(anti_var, .01)
        if var_ratio > 1:
            logging.info("Targets are %.2f x more variable than antitargets",
                         var_ratio)
        else:
            logging.info("Antitargets are %.2f x more variable than targets",
                         1. / var_ratio)

    if ((ref_matched[spread_key] > epsilon).any() and
        (np.abs(np.mod(ref_matched[log2_key], 1)) > epsilon).any()):
        # Pooled/paired reference only
        logging.debug("Weighting bins by coverage spread in reference")
        # NB: spread ~= SD, so variance ~= spread^2
        fancy_wt = 1.0 - ref_matched[spread_key] ** 2
        # Average w/ simple weights, giving this more emphasis
        x = .9
        weights = (x * fancy_wt
                   + (1 - x) * simple_wt)
    else:
        # Flat reference, only 1 weight estimate
        weights = simple_wt

    return cnarr.add_columns(weight=weights.clip(epsilon, 1.0))
