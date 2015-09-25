"""Supporting functions for the 'fix' command."""
from __future__ import absolute_import, division, print_function
import bisect

import numpy as np
import pandas as pd
from Bio._py3k import zip

from . import params, smoothing
from .ngfrills import echo


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
    echo("Keeping", sum(ok_cvg_indices), "of", len(ref_matched), "bins")
    pset = pset[ok_cvg_indices]
    ref_matched = ref_matched[ok_cvg_indices]

    # Apply corrections for known systematic biases in coverage
    pset.center_all()
    if fix_gc:
        if 'gc' in ref_matched:
            echo("Correcting for GC bias...")
            pset = center_by_window(pset, .1, ref_matched['gc'])
        else:
            echo("WARNING: Skipping correction for RepeatMasker bias")
    if fix_edge:
        echo("Correcting for density bias...")
        pset = center_by_window(pset, .1,
                                make_edge_sorter(pset, params.INSERT_SIZE))
    if fix_rmask:
        if 'rmask' in ref_matched:
            echo("Correcting for RepeatMasker bias...")
            pset = center_by_window(pset, .1, ref_matched['rmask'])
        else:
            echo("WARNING: Skipping correction for RepeatMasker bias")

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
    if isinstance(sort_key, (np.ndarray, pd.Series)):
        # Apply the same shuffling to the key array as to the target probe set
        sort_key = sort_key[shuffle_order]
    elif callable(sort_key):
        sort_key = df.apply(sort_key, axis=1) # ???
    else:
        raise ValueError("What is this?: %r" % sort_key)
    # Sort the data according to the specified parameter
    order = np.argsort(sort_key, kind='mergesort')
    df = df.iloc[order]
    biases = smoothing.rolling_median(df['log2'], fraction)
    # biases = smoothing.smoothed(df['log2'], fraction)
    df['log2'] -= biases
    fixarr = cnarr.as_dataframe(df)
    fixarr.sort()
    return fixarr


def make_edge_sorter(target_probes, margin):
    """Create a sort-key function for tiling edge effects."""
    # Index the target interval positions
    chrom_tile_starts = {}
    chrom_tile_ends = {}
    for chrom, rows in target_probes.by_chromosome():
        chrom_tile_starts[chrom] = rows['start']
        chrom_tile_ends[chrom] = rows['end']

    def get_edge(chrom, tgt_start, tgt_end, insert_size):
        """Quantify the "edge effect" of the target tile and its neighbors.

        The result is proportional to the change in the target's coverage due to
        these edge effects, i.e. the expected loss of coverage near the target
        edges and, if there are close neighboring tiles, gain of coverage due
        to "spill over" reads from the neighbor tiles.

        (This is not the actual change in coverage. This is just a tribute.)
        """
        margin_start = tgt_start - insert_size
        margin_end = tgt_end + insert_size
        tile_starts = chrom_tile_starts[chrom]
        tile_ends = chrom_tile_ends[chrom]
        target_size = (tgt_end - tgt_start)

        # Calculate coverage loss at (both) tile edges
        loss = edge_loss(target_size, insert_size)

        # For each neighbor tile, calculate coverage gain to the target
        gaps_left = []
        gaps_right = []
        # Find the leftmost tile in the margin
        left_idx = max(0, bisect.bisect_left(tile_ends, margin_start) - 1)
        for (tile_start, tile_end) in zip(tile_starts[left_idx:],
                                          tile_ends[left_idx:]):
            if tile_end <= margin_start:
                # No overlap on the 5' end -- keep moving forward
                continue
            if tile_start >= margin_end:
                # No overlap on the 3' end -- we're done
                break
            if tile_start == tgt_start and tile_end == tgt_end:
                # The target itself
                continue
            # Tile is within margins
            if margin_start <= tile_end <= tgt_start:
                # Left neighbor
                gaps_left.append(tgt_start - tile_end)
            elif tgt_end <= tile_start <= margin_end:
                # Right neighbor
                gaps_right.append(tile_start - tgt_end)
            elif tile_start < tgt_start and tile_end >= tgt_start:
                # Overlap on left side -- treat as adjacent
                gaps_left.append(0)
            elif tile_start <= tgt_end and tile_end > tgt_end:
                # Overlap on right side -- treat as adjacent
                gaps_right.append(0)
            else:
                # DBG: This should probably never happen
                echo("Oddly positioned tile (%s:%d-%d) vs. target (%d-%d)"
                     % (chrom, tile_start, tile_end, tgt_start, tgt_end))
                continue
        gain = 0
        if gaps_left:
            gain += edge_gain(target_size, insert_size, min(gaps_left))
        if gaps_right:
            gain += edge_gain(target_size, insert_size, min(gaps_right))
        return gain - loss

    def sorter_edge(row):
        """Calculate the edge effects on this bin.

        Find tiled intervals within a margin (+/- bp) of the given probe
        (including the probe itself, so the edge is never zero). Return the
        proportion of the windowed range that is covered by tiled regions.
        """
        return get_edge(row['chromosome'], row['start'], row['end'], margin)

    return sorter_edge


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
    """
    assert gap_size <= insert_size
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
        echo("Weighting bins by relative coverage depths in reference")
        # Penalize bins that deviate from expected coverage
        flat_cvgs = ref_matched.expect_flat_cvg()
        weights *= 2 ** -np.abs(ref_matched['log2'] - flat_cvgs)
    if (ref_matched['spread'] > epsilon).any():
        # NB: Not used with a flat or paired reference
        echo("Weighting bins by coverage spread in reference")
        # Inverse of variance, 0--1
        variances = ref_matched['spread'] ** 2
        invvars = 1.0 - (variances / variances.max())
        weights = (weights + invvars) / 2
    # Avoid 0-value bins -- CBS doesn't like these
    weights = np.maximum(weights, epsilon)
    return cnarr.add_columns(weight=weights)

