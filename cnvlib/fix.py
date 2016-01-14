"""Supporting functions for the 'fix' command."""
from __future__ import absolute_import, division, print_function

import logging

import numpy as np
import pandas as pd
from Bio._py3k import zip

from . import params, smoothing


def load_adjust_coverages(pset, ref_pset,
                          fix_gc, fix_edge, fix_rmask):
    """Load and filter probe coverages; correct using reference and GC."""
    if 'gc' in pset:
        # Don't choke on Picard-derived files that have the GC column
        pset = pset.drop_extra_columns()

    # No corrections needed if there are no data rows (e.g. no antitargets)
    if not len(pset):
        return pset

    ref_matched = match_ref_to_probes(ref_pset, pset)

    # Drop probes that had poor coverage in the pooled reference
    ok_cvg_indices = ~mask_bad_probes(ref_matched)
    logging.info("Keeping %d of %d bins", sum(ok_cvg_indices), len(ref_matched))
    pset = pset[ok_cvg_indices]
    ref_matched = ref_matched[ok_cvg_indices]

    # Apply corrections for known systematic biases in coverage
    pset.center_all()
    if fix_gc:
        if 'gc' in ref_matched:
            logging.info("Correcting for GC bias...")
            pset = center_by_window(pset, .1, ref_matched['gc'])
        else:
            logging.warn("WARNING: Skipping correction for RepeatMasker bias")
    if fix_edge:
        logging.info("Correcting for density bias...")
        edge_bias = get_edge_bias(pset, params.INSERT_SIZE)
        pset = center_by_window(pset, .1, edge_bias)
    if fix_rmask:
        if 'rmask' in ref_matched:
            logging.info("Correcting for RepeatMasker bias...")
            pset = center_by_window(pset, .1, ref_matched['rmask'])
        else:
            logging.warn("WARNING: Skipping correction for RepeatMasker bias")

    # Normalize coverages according to the reference
    # (Subtract the reference log2 copy number to get the log2 ratio)
    pset.data['log2'] -= ref_matched['log2']

    pset.center_all()
    return apply_weights(pset, ref_matched)


def mask_bad_probes(probes):
    """Flag the probes with excessively low or inconsistent coverage.

    Returns a bool array where True indicates probes that failed the checks.
    """
    mask = ((probes['log2'] < params.MIN_REF_COVERAGE) |
            (probes['spread'] > params.MAX_REF_SPREAD))
    if 'rmask' in probes:
        mask |= (probes['rmask'] > params.MAX_REPEAT_FRACTION)
    return mask


def match_ref_to_probes(ref_pset, probes):
    """Filter the reference probes to match the target or antitarget probe set.
    """
    probes_labeled = probes.data.set_index(probes.labels())
    ref_labeled = ref_pset.data.set_index(ref_pset.labels())
    # Safety
    for dset, name in ((probes_labeled, "probe"),
                       (ref_labeled, "reference")):
        dupes = dset.index.duplicated()
        if dupes.any():
            raise ValueError("Duplicated genomic coordinates in " + name +
                             " set:\n" + "\n".join(map(str, dset.index[dupes])))
    ref_matched = ref_labeled.reindex(index=probes_labeled.index)
    # Check for signs that the wrong reference was used
    num_missing = pd.isnull(ref_matched.start).sum()
    if num_missing > 0:
        raise ValueError("Reference is missing %d bins found in %s"
                         % (num_missing, probes.sample_id))
    return ref_pset.as_dataframe(ref_matched.reset_index(drop=True))


def center_by_window(cnarr, fraction, sort_key):
    """Smooth out biases according to the trait specified by sort_key.

    E.g. correct GC-biased probes by windowed averaging across similar-GC
    probes; or for similar interval sizes.
    """
    # Separate neighboring probes that could have the same key
    # (to avoid re-centering actual CNV regions -- only want an independently
    # sampled subset of presumably overall-CN-neutral probes)
    df = cnarr.data.reset_index(drop=True)
    shuffle_order = np.random.permutation(df.index)
    df = df.reindex(shuffle_order)
    # Apply the same shuffling to the key array as to the target probe set
    assert isinstance(sort_key, (np.ndarray, pd.Series))
    sort_key = sort_key[shuffle_order]
    # Sort the data according to the specified parameter
    order = np.argsort(sort_key, kind='mergesort')
    df = df.iloc[order]
    biases = smoothing.rolling_median(df['log2'], fraction)
    # biases = smoothing.smoothed(df['log2'], fraction)
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
        # tile_starts = subarr['start']
        # tile_ends = subarr['end']
        table = subarr.data
        table["margin_start"] = table.start - margin
        table["margin_end"] = table.end + margin
        table["target_size"] = table.end - table.start

        # Calculate coverage loss at (both edges of) each tile
        losses = table.target_size.apply(edge_loss, args=(margin,))

        # Find the neighboring tiles within each tile's margins,
        # but not the target itself
        table["l_start_idx"] = table.end.searchsorted(table.margin_start)
        table["l_end_idx"] = table.start.searchsorted(table.start)
        table["r_start_idx"] = table.end.searchsorted(table.end, "right")
        table["r_end_idx"] = table.start.searchsorted(table.margin_end, "right")

        # For each target, calculate coverage gain from neighboring targets
        # Find tiled intervals within a margin (+/- bp) of the given probe
        # (excluding the probe itself), then calculate the relative coverage
        # "gain" due to the neighbors, if any.
        def l_row_gains(row):
            gaps_left = row.start - table.end[row.l_start_idx:row.l_end_idx]
            return edge_gain(row.target_size, margin, gaps_left.min())

        def r_row_gains(row):
            gaps_right = table.start[row.r_start_idx:row.r_end_idx] - row.end
            return edge_gain(row.target_size, margin, gaps_right.min())

        gains = np.zeros(len(table))
        l_rows = table[table.l_start_idx < table.l_end_idx]
        gains[l_rows.index.values] = l_rows.apply(l_row_gains,
                                                  axis=1, reduce=True)

        r_rows = table[table.r_start_idx < table.r_end_idx]
        gains[r_rows.index.values] += r_rows.apply(r_row_gains,
                                                   axis=1, reduce=True)

        out_row = gains - losses
        output_by_chrom.append(out_row)
    return np.concatenate(output_by_chrom)


def edge_loss(target_size, insert_size):
    """Calculate coverage loss at the edges of a baited region.

    Letting i = insert size and t = target size, the proportional loss of
    coverage near the two edges of the baited region (combined) is::

        i/2t

    If the "shoulders" extend outside the bait $(t < i), reduce by::

        (i-t)^2 / 4it

    on each side, or (i-t)^2 / 2it total.
    """
    loss = insert_size / (2 * target_size)
    if target_size < insert_size:
        # Drop the shoulder part that would extend past the bait
        loss -= ((insert_size - target_size)**2
                 / (2 * insert_size * target_size))
    return loss


def edge_gain(target_size, insert_size, gap_size):
    """Calculate coverage gain from a neighboring bait's flanking reads.

    Letting i = insert size, t = target size, g = gap to neighboring bait,
    the gain of coverage due to a nearby bait, if g < i, is::

        (i-g)^2 / 4it

    If the neighbor flank extends beyond the target (t+g < i), reduce by::

        (i-t-g)^2 / 4it

    If a neighbor overlaps the target, treat it as adjacent (gap size 0).
    """
    assert gap_size <= insert_size, (gap_size, insert_size)
    gap_size = max(0, gap_size)
    gain = ((insert_size - gap_size)**2
            / (4 * insert_size * target_size))
    if target_size + gap_size < insert_size:
        # Drop the flank part that extends past this baited region
        gain -= ((insert_size - target_size - gap_size)**2
                 / (4 * insert_size * target_size))
    return gain


def apply_weights(cnarr, ref_matched, epsilon=1e-4):
    """Calculate weights for each bin.

    Weights are derived from:

    - bin sizes
    - average bin coverage depths in the reference
    - the "spread" column of the reference.
    """
    # Relative bin sizes
    sizes = ref_matched['end'] - ref_matched['start']
    weights = sizes / sizes.max()
    if (np.abs(np.mod(ref_matched['log2'], 1)) > epsilon).any():
        # NB: Not used with a flat reference
        logging.info("Weighting bins by relative coverage depths in reference")
        # Penalize bins that deviate from expected coverage
        flat_cvgs = ref_matched.expect_flat_cvg()
        weights *= 2 ** -np.abs(ref_matched['log2'] - flat_cvgs)
    if (ref_matched['spread'] > epsilon).any():
        # NB: Not used with a flat or paired reference
        logging.info("Weighting bins by coverage spread in reference")
        # Inverse of variance, 0--1
        variances = ref_matched['spread'] ** 2
        invvars = 1.0 - (variances / variances.max())
        weights = (weights + invvars) / 2
    # Avoid 0-value bins -- CBS doesn't like these
    weights = np.maximum(weights, epsilon)
    return cnarr.add_columns(weight=weights)

