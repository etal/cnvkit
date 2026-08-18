"""Plotting utilities."""

from __future__ import annotations

import collections
import logging
import math
from typing import TYPE_CHECKING

import numpy as np

from skgenome.intersect import require_ascending_starts
from skgenome.rangelabel import Region, unpack_range

from . import cnary, params

if TYPE_CHECKING:
    from matplotlib.axes._axes import Axes

    from cnvlib.cnary import CopyNumArray


MB = 1e-6  # To rescale from bases to megabases


def binwise_x(varr):
    """X positions of `varr` on the bin axis, sub-bin offsets included.

    `update_binwise_positions` spreads variants that share a bin across its
    width, and parks that offset in a ``bin_offset`` column rather than in
    `start`: it is a rendering artifact with no genomic meaning, and
    `GenomicArray` recasts its coordinate columns to int on every rewrap, so a
    fraction stored in `start` would not survive to be drawn. A caller's own
    ``bin_offset`` column, were one ever read from a file, is overwritten.

    Arrays that never went through that transform carry no such column and
    plot at `start`, so callers need not know which axis they are on.
    """
    x = varr["start"].values
    if "bin_offset" in varr:
        return x + varr["bin_offset"].values
    return x


def plot_chromosome_dividers(
    axis: Axes,
    chrom_sizes: collections.OrderedDict,
    pad: float | None = None,
    along: str = "x",
) -> collections.OrderedDict:
    """Given chromosome sizes, plot divider lines and labels.

    Draws black lines between each chromosome, with padding. Labels each chromosome range with the chromosome name,
    centered in the region, under a tick. Sets the axis limits to the covered range.

    By default, the dividers are vertical and the labels are on the X axis of the plot. If the `along` parameter is 'y',
    this is transposed to horizontal dividers and the labels on the Y axis.

    Returns
    -------
    OrderedDict
        A table of the position offsets of each chromosome along the specified axis.
    """
    assert isinstance(chrom_sizes, collections.OrderedDict)
    if pad is None:
        pad = 0.003 * sum(chrom_sizes.values())
    dividers = []
    centers = []
    starts = collections.OrderedDict()
    curr_offset = pad
    for label, size in list(chrom_sizes.items()):
        starts[label] = curr_offset
        centers.append(curr_offset + 0.5 * size)
        dividers.append(curr_offset + size + pad)
        curr_offset += size + 2 * pad

    if along not in ("x", "y"):
        raise ValueError(
            "Direction for plotting chromosome dividers and labels along must be either x or y."
        )

    if along == "x":
        axis.set_xlim(0, curr_offset)
        for xposn in dividers[:-1]:
            axis.axvline(x=xposn, color="k")
        # Use chromosome names as x-axis labels (instead of base positions)
        axis.set_xticks(centers)
        axis.set_xticklabels(list(chrom_sizes.keys()), rotation=90)
        axis.tick_params(labelsize="small")
        axis.tick_params(axis="x", length=0)
        axis.get_yaxis().tick_left()
    else:
        axis.set_ylim(0, curr_offset)
        for yposn in dividers[:-1]:
            axis.axhline(y=yposn, color="k")
        # Use chromosome names as y-axis labels (instead of base positions)
        axis.set_yticks(centers)
        axis.set_yticklabels(list(chrom_sizes.keys()))
        axis.tick_params(labelsize="small")
        axis.tick_params(axis="y", length=0)
        axis.get_xaxis().tick_bottom()

    return starts


# ________________________________________
# Internal supporting functions


def chromosome_sizes(
    probes: CopyNumArray, to_mb: bool = False
) -> collections.OrderedDict:
    """Create an ordered mapping of chromosome names to sizes."""
    chrom_sizes = collections.OrderedDict()
    for chrom, rows in probes.by_chromosome():
        # Use the rightmost coordinate from either column, so an array holding
        # a reverse-oriented interval (start > end) still sets the
        # chromosome's extent. `skgenome.tabio.read` keeps such a row out of
        # every region and copy-number array it returns, but arrays built by
        # library callers are not covered. For well-formed intervals end.max()
        # already dominates, so this is a no-op there.
        chrom_sizes[chrom] = max(rows["end"].max(), rows["start"].max())
        if to_mb:
            chrom_sizes[chrom] *= MB
    return chrom_sizes


def translate_region_to_bins(region, bins):
    """Map genomic coordinates to bin indices.

    Return a tuple of (chrom, start, end), just like unpack_range.

    `bins` must hold the named chromosome's rows in ascending start order: a
    bin index is a position in that order, so a coordinate's rank among the
    bin starts is only its index on the axis while the two agree.
    """
    if region is None:
        return Region(None, None, None)
    chrom, start, end = unpack_range(region)
    if start is None and end is None:
        return Region(chrom, start, end)
    if start is None:
        start = 0
    if end is None:
        end = float("inf")
    c_rows = bins.data.loc[bins.data.chromosome == chrom]
    require_ascending_starts(c_rows)
    # NB: only bin start positions matter here
    r_start, r_end = np.searchsorted(c_rows["start"].values, [start, end])
    return Region(chrom, r_start, r_end)


def update_binwise_positions_simple(cnarr):
    start_chunks = []
    end_chunks = []
    is_segment = "probes" in cnarr
    if is_segment:
        cnarr = cnarr[cnarr["probes"] > 0]
    for _chrom, c_arr in cnarr.by_chromosome():
        if is_segment:
            # Segments -- each row can cover many bins
            ends = c_arr["probes"].values.cumsum()
            starts = np.r_[0, ends[:-1]]
        else:
            # Bins -- enumerate rows
            n_bins = len(c_arr)
            starts = np.arange(n_bins)
            ends = np.arange(1, n_bins + 1)
        start_chunks.append(starts)
        end_chunks.append(ends)
    return cnarr.as_dataframe(
        cnarr.data.assign(
            start=np.concatenate(start_chunks), end=np.concatenate(end_chunks)
        )
    )


def update_binwise_positions(
    cnarr, segments=None, variants=None, *, extra_variants=None
):
    """Convert start/end positions from genomic to bin-wise coordinates.

    Instead of chromosomal basepairs, the positions indicate enumerated bins.

    The bin axis is `cnarr`'s row order -- `np.arange` over each chromosome's
    rows -- while `searchsorted` below returns a coordinate's rank among the
    bin starts. Those coincide only while each chromosome's rows ascend, and
    there is no answer to compute when they do not: a correct rank drawn
    against a scrambled axis is still a wrong figure. So rows out of order
    raise `ValueError` rather than rendering silently.

    Revise the start and end values for all GenomicArray instances at once,
    where the `cnarr` bins are mapped to corresponding `segments`, and
    `variants` are placed on the same axis.

    A variant inside a bin takes that bin's ordinal. One between two bins has
    no bin of its own, and the gap has no width on this axis, so it sits on
    the boundary. One before the first bin or past the last has nowhere to
    sit and is dropped: on the genomic axis its distance from any covered bin
    is visible, here it would sit on a real bin and claim a place in the data
    it never had.

    `start` stays an integer ordinal; variants sharing a bin are separated by
    an offset that `binwise_x` adds back at draw time.

    Each variant also gains a ``genomic_segment`` column naming the segment
    that contains it, resolved here because this is the last moment the
    genomic coordinates exist; see `_freeze_segment_membership`.

    ``extra_variants`` is an optional list of additional VariantArray (or None)
    instances to translate against the same bin enumeration as ``variants``
    -- used by the scatter plot's ``--show-snvs`` overlays (#290) so the LOH
    and somatic markers stay aligned with the het dots under ``--by-bin``.

    Returns a 4-tuple ``(cnarr, segments, variants, extras)`` where ``extras``
    is the list of translated extra arrays in input order (empty list if no
    extras were passed). ``variants`` and each of ``extras`` may hold fewer
    rows than the array they came from, the ones with no place on the bin
    axis having been dropped as described above.
    """
    cnarr = cnarr.copy()
    if segments:
        segments = segments.copy()
        seg_chroms = set(segments.chromosome.unique())
    # Bundle the primary `variants` and any extras so the per-chrom loop
    # iterates them uniformly. Copy each so the caller's arrays are unmodified.
    var_arrs = [variants, *(list(extra_variants) if extra_variants else [])]
    var_arrs = [v.copy() if v is not None else None for v in var_arrs]
    # Only the arrays that exist take part in the loop, each carrying the
    # chromosomes it spans and a mask of the rows that found a place on the
    # axis; the rest are dropped at the end. Their slots in `var_arrs` are
    # kept so the return preserves the caller's argument order.
    present = [
        (i, v, set(v.chromosome.unique()), np.zeros(len(v), dtype=bool))
        for i, v in enumerate(var_arrs)
        if v is not None
    ]
    if segments:
        for _i, varr, _c, _k in present:
            _freeze_segment_membership(varr, segments)

    # ENH: look into pandas groupby innards to get group indices
    for chrom in cnarr.chromosome.unique():
        # Enumerate bins, starting from 0
        # NB: plotted points will be at +0.5 offsets
        c_idx = cnarr.chromosome == chrom
        c_bins = cnarr[c_idx]  # .copy()
        require_ascending_starts(c_bins.data)
        bin_starts = c_bins.start.values
        bin_ends = c_bins.end.values
        if segments and chrom in seg_chroms:
            # Match segment boundaries to enumerated bins. Read each segment's
            # own end rather than chaining to the next segment's start, which
            # is wrong when the segments do not ascend and tile the bins.
            c_seg_idx = (segments.chromosome == chrom).values
            seg_starts = np.searchsorted(bin_starts, segments.start.values[c_seg_idx])
            seg_ends = np.searchsorted(bin_starts, segments.end.values[c_seg_idx])
            segments.data.loc[c_seg_idx, "start"] = seg_starts
            segments.data.loc[c_seg_idx, "end"] = seg_ends

        for _i, varr, vchroms, keep in present:
            if chrom not in vchroms:
                continue
            c_varr_idx = np.flatnonzero((varr.chromosome == chrom).values)
            positions = varr.start.values[c_varr_idx]
            # searchsorted gives the rank among bin starts, which is one past
            # the bin holding the position unless it is a bin start exactly.
            below = np.searchsorted(bin_starts, positions, side="right") - 1
            after_first_bin = below >= 0
            inside = after_first_bin & (positions < bin_ends[np.maximum(below, 0)])
            # In a gap: no bin of its own, so it sits on the boundary with the
            # next one. Past the last bin that boundary is off the axis.
            v_starts = np.where(inside, below, below + 1)
            on_axis = after_first_bin & (v_starts < len(bin_starts))

            rows = c_varr_idx[on_axis]
            keep[rows] = True
            v_starts = v_starts[on_axis]
            labels = varr.data.index[rows]
            varr.data.loc[labels, "start"] = v_starts
            # One bin wide, matching how the bins themselves are numbered.
            varr.data.loc[labels, "end"] = v_starts + 1
            varr.data.loc[labels, "bin_offset"] = _within_bin_offsets(v_starts)

        c_starts = np.arange(len(c_bins))  # c_idx.sum())
        c_ends = np.arange(1, len(c_bins) + 1)
        cnarr.data.loc[c_idx, "start"] = c_starts
        cnarr.data.loc[c_idx, "end"] = c_ends

    dropped = sum(int((~keep).sum()) for _i, _v, _c, keep in present)
    if dropped:
        total = sum(len(keep) for _i, _v, _c, keep in present)
        logging.info(
            "Dropped %d of %d variant(s) lying outside the binned regions; "
            "the bin axis has no position for them",
            dropped,
            total,
        )
    for i, varr, _c, keep in present:
        var_arrs[i] = varr[keep]
    return cnarr, segments, var_arrs[0], var_arrs[1:]


def _freeze_segment_membership(varr, segments) -> None:
    """Freeze which segment genomically contains each variant, in `varr`.

    The bin ordinal cannot answer this later. A variant in a gap between bins
    is parked on the boundary, which can lie inside a segment its genomic
    position never reached -- lending that segment a B-allele frequency it has
    no evidence for.

    The answer is the segment's index label, not its row number: `in_range`
    hands the drawing code a trimmed subset whose rows keep their labels but
    no longer sit at their original positions. Variants in no segment keep the
    -1 fill, which matches no label and so joins nothing. A caller's own
    ``genomic_segment`` column, were one ever read from a file, is overwritten.

    One label per variant, so a multi-base variant overlapping two segments
    joins the later one rather than both, where the genomic axis counts it
    twice. Segments CNVkit produces are point-disjoint, and no fixture has
    such a variant; the two axes otherwise agree.
    """
    labels = segments.data.index
    if not labels.is_unique or (labels == -1).any():
        raise ValueError(
            "Segments must carry unique index labels, none of them -1, to key "
            f"B-allele frequency membership: {labels.tolist()}"
        )
    varr.data["genomic_segment"] = -1
    for label, (_seg, group) in zip(labels, varr.by_ranges(segments), strict=True):
        if len(group):
            varr.data.loc[group.data.index, "genomic_segment"] = label


def _within_bin_offsets(bin_indices):
    """Spread variants that share a bin evenly across that bin's width.

    Returns one value in [0, 1) per entry: (2j+1)/2n for the j-th of n sharing
    a bin, so a variant alone in its bin lands at 0.5 -- the bin midpoint,
    where that bin's own log2 dot is drawn on the panel above.

    Grouped by value rather than by consecutive run, so the point cloud does
    not depend on the order the rows arrive in; a stable sort settles which
    row gets which slot.
    """
    order = np.argsort(bin_indices, kind="stable")
    ranked = bin_indices[order]
    starts_run = np.flatnonzero(np.concatenate(([True], ranked[1:] != ranked[:-1])))
    sizes = np.diff(np.append(starts_run, len(ranked)))
    offsets = np.empty(len(ranked))
    offsets[order] = (
        np.arange(len(ranked)) - np.repeat(starts_run, sizes) + 0.5
    ) / np.repeat(sizes, sizes)
    return offsets


# ________________________________________
# Utilies used by other modules


def cvg2rgb(cvg: float, desaturate: bool) -> tuple[float, float, float]:
    """Choose a shade of red or blue representing log2-coverage value."""
    cutoff = 1.33  # Values above this magnitude are shown with max intensity
    x = min(abs(cvg) / cutoff, 1.0)
    if desaturate:
        # Adjust intensity sigmoidally -- reduce near 0, boost near 1
        # Exponent <1 shifts the fixed point leftward (from x=0.5)
        x = ((1.0 - math.cos(x * math.pi)) / 2.0) ** 0.8
        # Slight desaturation of colors at lower coverage
        s = x**1.2
    else:
        s = x
    if cvg < 0:
        rgb = (1 - s, 1 - s, 1 - 0.25 * x)  # Blueish
    else:
        rgb = (1 - 0.25 * x, 1 - s, 1 - s)  # Reddish
    return rgb


def gene_coords_by_name(probes, names):
    """Find the chromosomal position of each named gene in probes.

    Each locus of a requested name gets its own region, grouped exactly as
    `CopyNumArray.by_gene` groups bins, via `cnvlib.cnary.gene_runs`. A name
    that recurs at several loci -- a repeat family, or a gene name reused on
    another chromosome -- therefore yields one region per locus instead of one
    span reaching from its first occurrence to its last.

    Returns
    -------
    dict
        Of: {chromosome: [(start, end, gene name), ...]}
    """
    names = set(filter(None, names))
    if not names:
        return {}

    ignore = (*params.IGNORE_GENE_NAMES, *params.ANTITARGET_ALIASES)
    # A locus spans its run's first and last bin -- the same start and end
    # 'genemetrics' publishes for that locus (reports.group_by_genes).
    # Regions are keyed so that requested names sharing a locus, e.g. the
    # co-binned pair "ERBB2,MIR4728", are labeled once; co-binned names the
    # caller did not ask for are never surfaced (#458).
    all_coords: dict = collections.defaultdict(lambda: collections.defaultdict(set))
    found: set[str] = set()
    for chrom, subarr in probes.by_chromosome():
        starts = subarr["start"].to_numpy()
        ends = subarr["end"].to_numpy()
        for start_pos, end_pos, gene in cnary.gene_runs(subarr["gene"], ignore):
            if gene in names:
                found.add(gene)
                all_coords[chrom][starts[start_pos], ends[end_pos - 1]].add(gene)
    if unfound := names - found:
        raise ValueError(
            "No targeted gene named " + ", ".join(map(repr, sorted(unfound))) + " found"
        )
    # Consolidate each region's gene names into a string
    uniq_coords = {}
    for chrom, hits in all_coords.items():
        uniq_coords[chrom] = [
            (start, end, ",".join(sorted(gene_names)))
            for (start, end), gene_names in hits.items()
        ]
    return uniq_coords


def gene_coords_by_range(probes, chrom, start, end, ignore=params.IGNORE_GENE_NAMES):
    """Find the chromosomal position of all genes in a range.

    Returns
    -------
    dict
        Of: {chromosome: [(start, end, gene), ...]}
    """
    ignore += params.ANTITARGET_ALIASES
    # Tabulate the genes in the selected region
    genes: dict = collections.OrderedDict()
    for row in probes.in_range(chrom, start, end):
        if not isinstance(row.gene, str):
            continue
        name = row.gene
        if name in genes:
            genes[name][1] = row.end
        elif name not in ignore:
            genes[name] = [row.start, row.end]
    # Reorganize the data structure
    return {
        chrom: [(gstart, gend, name) for name, (gstart, gend) in list(genes.items())]
    }
