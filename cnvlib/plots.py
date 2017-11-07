"""Plotting utilities."""
from __future__ import absolute_import, division, print_function
from builtins import str

import collections
import itertools
import logging
import math

import numpy as np

from . import core, params
from skgenome.rangelabel import unpack_range, Region


MB = 1e-6  # To rescale from bases to megabases


def plot_x_dividers(axis, chrom_sizes, pad=None):
    """Plot vertical dividers and x-axis labels given the chromosome sizes.

    Draws vertical black lines between each chromosome, with padding.
    Labels each chromosome range with the chromosome name, centered in the
    region, under a tick.
    Sets the x-axis limits to the covered range.

    Returns
    -------
    OrderedDict
        A table of the x-position offsets of each chromosome.
    """
    assert isinstance(chrom_sizes, collections.OrderedDict)
    if pad is None:
        pad = 0.003 * sum(chrom_sizes.values())
    x_dividers = []
    x_centers = []
    x_starts = collections.OrderedDict()
    curr_offset = pad
    for label, size in list(chrom_sizes.items()):
        x_starts[label] = curr_offset
        x_centers.append(curr_offset + 0.5 * size)
        x_dividers.append(curr_offset + size + pad)
        curr_offset += size + 2 * pad

    axis.set_xlim(0, curr_offset)
    for xposn in x_dividers[:-1]:
        axis.axvline(x=xposn, color='k')
    # Use chromosome names as x-axis labels (instead of base positions)
    axis.set_xticks(x_centers)
    axis.set_xticklabels(list(chrom_sizes.keys()), rotation=60)
    axis.tick_params(labelsize='small')
    axis.tick_params(axis='x', length=0)
    axis.get_yaxis().tick_left()

    return x_starts


# ________________________________________
# Internal supporting functions

def chromosome_sizes(probes, to_mb=False):
    """Create an ordered mapping of chromosome names to sizes."""
    chrom_sizes = collections.OrderedDict()
    for chrom, rows in probes.by_chromosome():
        chrom_sizes[chrom] = rows['end'].max()
        if to_mb:
            chrom_sizes[chrom] *= MB
    return chrom_sizes


def translate_region_to_bins(region, bins):
    """Map genomic coordinates to bin indices.

    Return a tuple of (chrom, start, end), just like unpack_range.
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
    # NB: only bin start positions matter here
    c_bin_starts = bins.data.loc[bins.data.chromosome == chrom, "start"].values
    r_start, r_end = np.searchsorted(c_bin_starts, [start, end])
    return Region(chrom, r_start, r_end)


def translate_segments_to_bins(segments, bins):
    if "probes" in segments and segments["probes"].sum() == len(bins):
        # Segments and .cnr bins already match
        return update_binwise_positions_simple(segments)
    else:
        logging.warning("Segments %s 'probes' sum does not match the number "
                        "of bins in %s", segments.sample_id, bins.sample_id)
        # Must re-align segments to .cnr bins
        _x, segments, _v = update_binwise_positions(bins, segments)
        return segments


def update_binwise_positions_simple(cnarr):
    start_chunks = []
    end_chunks = []
    is_segment = ("probes" in cnarr)
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
        cnarr.data.assign(start=np.concatenate(start_chunks),
                          end=np.concatenate(end_chunks)))


def update_binwise_positions(cnarr, segments=None, variants=None):
    """Convert start/end positions from genomic to bin-wise coordinates.

    Instead of chromosomal basepairs, the positions indicate enumerated bins.

    Revise the start and end values for all GenomicArray instances at once,
    where the `cnarr` bins are mapped to corresponding `segments`, and
    `variants` are grouped into `cnarr` bins as well -- if multiple `variants`
    rows fall within a single bin, equally-spaced fractional positions are used.

    Returns copies of the 3 input objects with revised `start` and `end` arrays.
    """
    cnarr = cnarr.copy()
    if segments:
        segments = segments.copy()
        seg_chroms = set(segments.chromosome.unique())
    if variants:
        variants = variants.copy()
        var_chroms = set(variants.chromosome.unique())

    # ENH: look into pandas groupby innards to get group indices
    for chrom in cnarr.chromosome.unique():
        # Enumerate bins, starting from 0
        # NB: plotted points will be at +0.5 offsets
        c_idx = (cnarr.chromosome == chrom)
        c_bins = cnarr[c_idx]#.copy()
        if segments and chrom in seg_chroms:
            # Match segment boundaries to enumerated bins
            c_seg_idx = (segments.chromosome == chrom).values
            seg_starts = np.searchsorted(c_bins.start.values,
                                         segments.start.values[c_seg_idx])
            seg_ends = np.r_[seg_starts[1:], len(c_bins)]
            segments.data.loc[c_seg_idx, "start"] = seg_starts
            segments.data.loc[c_seg_idx, "end"] = seg_ends

        if variants and chrom in var_chroms:
            # Match variant positions to enumerated bins, and
            # add fractional increments to multiple variants within 1 bin
            c_varr_idx = (variants.chromosome == chrom).values
            c_varr_df = variants.data[c_varr_idx]
            # Get binwise start indices of the variants
            v_starts = np.searchsorted(c_bins.start.values,
                                       c_varr_df.start.values)
            # Overwrite runs of repeats with fractional increments,
            #   adding the cumulative fraction to each repeat
            for idx, size in list(get_repeat_slices(v_starts)):
                v_starts[idx] += np.arange(size) / size
            variant_sizes = c_varr_df.end - c_varr_df.start
            variants.data.loc[c_varr_idx, "start"] = v_starts
            variants.data.loc[c_varr_idx, "end"] = v_starts + variant_sizes

        c_starts = np.arange(len(c_bins)) # c_idx.sum())
        c_ends = np.arange(1, len(c_bins) + 1)
        cnarr.data.loc[c_idx, "start"] = c_starts
        cnarr.data.loc[c_idx, "end"] = c_ends

    return cnarr, segments, variants


def get_repeat_slices(values):
    """Find the location and size of each repeat in `values`."""
    # ENH: look into pandas groupby innards
    offset = 0
    for idx, (_val, rpt) in enumerate(itertools.groupby(values)):
        size = len(list(rpt))
        if size > 1:
            i = idx + offset
            slc = slice(i, i+size)
            yield slc, size
            offset += size - 1



# ________________________________________
# Utilies used by other modules

def cvg2rgb(cvg, desaturate):
    """Choose a shade of red or blue representing log2-coverage value."""
    cutoff = 1.33  # Values above this magnitude are shown with max intensity
    x = min(abs(cvg) / cutoff, 1.0)
    if desaturate:
        # Adjust intensity sigmoidally -- reduce near 0, boost near 1
        # Exponent <1 shifts the fixed point leftward (from x=0.5)
        x = ((1. - math.cos(x * math.pi)) / 2.) ** 0.8
        # Slight desaturation of colors at lower coverage
        s = x**1.2
    else:
        s = x
    if cvg < 0:
        rgb = (1 - s, 1 - s, 1 - .25*x)  # Blueish
    else:
        rgb = (1 - .25*x, 1 - s, 1 - s)  # Reddish
    return rgb


# XXX should this be a CopyNumArray method?
# or: use by_genes internally
# or: have by_genes use this internally
def gene_coords_by_name(probes, names):
    """Find the chromosomal position of each named gene in probes.

    Returns
    -------
    dict
        Of: {chromosome: [(start, end, gene name), ...]}
    """
    names = list(filter(None, set(names)))
    if not names:
        return {}

    # Create an index of gene names
    gene_index = collections.defaultdict(set)
    for i, gene in enumerate(probes['gene']):
        for gene_name in gene.split(','):
            if gene_name in names:
                gene_index[gene_name].add(i)
    # Retrieve coordinates by name
    all_coords = collections.defaultdict(lambda : collections.defaultdict(set))
    for name in names:
        gene_probes = probes.data.take(sorted(gene_index.get(name, [])))
        if not len(gene_probes):
            raise ValueError("No targeted gene named '%s' found" % name)
        # Find the genomic range of this gene's probes
        start = gene_probes['start'].min()
        end = gene_probes['end'].max()
        chrom = core.check_unique(gene_probes['chromosome'], name)
        # Deduce the unique set of gene names for this region
        uniq_names = set()
        for oname in set(gene_probes['gene']):
            uniq_names.update(oname.split(','))
        all_coords[chrom][start, end].update(uniq_names)
    # Consolidate each region's gene names into a string
    uniq_coords = {}
    for chrom, hits in all_coords.items():
        uniq_coords[chrom] = [(start, end, ",".join(sorted(gene_names)))
                              for (start, end), gene_names in hits.items()]
    return uniq_coords


def gene_coords_by_range(probes, chrom, start, end,
                         ignore=params.IGNORE_GENE_NAMES):
    """Find the chromosomal position of all genes in a range.

    Returns
    -------
    dict
        Of: {chromosome: [(start, end, gene), ...]}
    """
    ignore += params.ANTITARGET_ALIASES
    # Tabulate the genes in the selected region
    genes = collections.OrderedDict()
    for row in probes.in_range(chrom, start, end):
        name = str(row.gene)
        if name in genes:
            genes[name][1] = row.end
        elif name not in ignore:
            genes[name] = [row.start, row.end]
    # Reorganize the data structure
    return {chrom: [(gstart, gend, name)
                    for name, (gstart, gend) in list(genes.items())]}
