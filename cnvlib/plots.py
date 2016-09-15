"""Plotting utilities."""
from __future__ import absolute_import, division
from builtins import str, zip
from past.builtins import basestring

import collections
import logging
import math
import sys

import numpy as np

from . import core, params, smoothing


SEG_COLOR = "darkorange"
POINT_COLOR = '#606060'

# === Chromosome-level scatter plots ===

HIGHLIGHT_COLOR = 'gold'
MB = 1e-6  # To rescale from bases to megabases

def setup_chromosome(axis, probes=None, segments=None, variants=None,
                     y_min=None, y_max=None, y_label=None):
    """Configure axes for plotting a single chromosome's data.

    Parameters
    ----------
    probes : CopyNumArray
    segments : CopyNumArray
    variants : VariantArray
        All should already be subsetted to the region that will be plotted.
    """
    min_x = np.inf
    max_x = 0
    for arr in (probes, segments, variants):
        if arr and len(arr):
            max_x = max(max_x, arr.end.iat[-1])
            min_x = min(min_x, arr.start.iat[0])
    if not max_x > min_x:
        if any((probes, segments, variants)):
            logging.warn("*WARNING* selection start %s > end %s (%s); did you "
                         "correctly sort the input file by genomic location?",
                         min_x, max_x, max_x - min_x)
        raise ValueError("No usable data points to plot out of "
                         "%d probes, %d segments, %d variants"
                         % (len(probes) if probes else 0,
                            len(segments) if segments else 0,
                            len(variants) if variants else 0))
    axis.set_xlim(min_x * MB, max_x * MB)
    if y_min and y_max:
        axis.set_ylim(y_min, y_max)
        if y_min < 0 < y_max:
            axis.axhline(color='k')
    if y_label:
        axis.set_ylabel(y_label)
    axis.tick_params(which='both', direction='out')
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()


def cnv_on_chromosome(axis, probes, segments, genes, background_marker=None,
                      do_trend=False, y_min=None, y_max=None):
    """Draw a scatter plot of probe values with CBS calls overlaid.

    Parameters
    ----------
    genes : list
        Of tuples: (start, end, gene name)
    """
    # Get scatter plot coordinates
    x = 0.5 * (probes['start'] + probes['end']) * MB # bin midpoints
    y = probes['log2']
    if 'weight' in probes:
        w = 46 * probes['weight'] ** 2 + 2
    else:
        w = np.repeat(30, len(x))
    is_bg = (probes['gene'] == 'Background')

    # Configure axes
    # TODO - use segment y-values if probes not given
    if not y_min:
        y_min = max(-5.0, min(y.min() - .1, -.3)) if len(y) else -1.1
    if not y_max:
        y_max = max(.3, y.max() + (.25 if genes else .1)) if len(y) else 1.1
    setup_chromosome(axis, probes, segments, None, y_min, y_max,
                     "Copy ratio (log2)")
    if genes:
        # Rotate text in proportion to gene density
        ngenes = len(genes)
        text_size = ('small' if ngenes <= 6 else 'x-small')
        if ngenes <= 3:
            text_rot = 'horizontal'
        elif ngenes <= 6:
            text_rot = 30
        elif ngenes <= 10:
            text_rot = 45
        elif ngenes <= 20:
            text_rot = 60
        else:
            text_rot = 'vertical'
        for gene in genes:
            gene_start, gene_end, gene_name = gene
            # Highlight and label gene region
            # (rescale positions from bases to megabases)
            axis.axvspan(gene_start * MB, gene_end * MB,
                         alpha=0.5, color=HIGHLIGHT_COLOR, zorder=-1)
            axis.text(0.5 * (gene_start + gene_end) * MB,
                      min(2.4, y.max() + .1) if len(y) else .1,
                      gene_name,
                      horizontalalignment='center',
                      rotation=text_rot,
                      size=text_size)
                      # size='small')

    if background_marker in (None, 'o'):
        # Plot targets and antitargets with the same marker
        axis.scatter(x, y, w, color=POINT_COLOR, alpha=0.4, marker='o')
    else:
        # Use the given marker to plot antitargets separately
        x_fg = []
        y_fg = []
        w_fg = []
        x_bg = []
        y_bg = []
        # w_bg = []
        for x_pt, y_pt, w_pt, is_bg_pt in zip(x, y, w, is_bg):
            if is_bg_pt:
                x_bg.append(x_pt)
                y_bg.append(y_pt)
                # w_bg.append(w_pt)
            else:
                x_fg.append(x_pt)
                y_fg.append(y_pt)
                w_fg.append(w_pt)
        axis.scatter(x_fg, y_fg, w_fg, color=POINT_COLOR, alpha=0.4, marker='o')
        axis.scatter(x_bg, y_bg, color=POINT_COLOR, alpha=0.5,
                     marker=background_marker)

    # Add a local trend line
    if do_trend:
        axis.plot(x, smoothing.smoothed(y, 50),
                    color=POINT_COLOR, linewidth=2, zorder=-1)

    # Get coordinates for CBS lines & draw them
    if segments:
        for row in segments:
            axis.plot((row.start * MB, row.end * MB),
                      (row.log2, row.log2),
                      color=SEG_COLOR, linewidth=4, solid_capstyle='round')


def snv_on_chromosome(axis, variants, segments, genes,
                      do_trend, do_boost=False):
    # XXX only set x-limits if not already done for probes/segments
    # setup_chromosome(axis, None, segments, variants,
    #                  0.0, 1.0, "VAF")
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("VAF")
    axis.set_xlabel("Position (Mb)")
    axis.get_yaxis().tick_left()
    axis.get_xaxis().tick_top()
    axis.tick_params(which='both', direction='out',
                     labelbottom=False, labeltop=False)

    x_mb = variants["start"] * MB
    if do_boost:
        y = variants.tumor_boost()
    else:
        y = np.asfarray(variants["alt_freq"])
    axis.scatter(x_mb, y, color=POINT_COLOR, alpha=0.3)
    # TODO - highlight genes/selection
    if segments:
        # Draw average VAF within each segment
        posns = np.asfarray(variants["start"]) # * MB
        y = np.asfarray(y)
        for v_start, v_end, v_freq in group_snvs_by_segments(posns, y,
                                                             segments):
            # ENH: color by segment gain/loss
            axis.plot([v_start * MB, v_end * MB], [v_freq, v_freq],
                      color='#C0C0C0', linewidth=2, #zorder=1,
                      solid_capstyle='round')

# === Genome-level scatter plots ===

def setup_genome(axis, probes, segments, variants, y_min=None, y_max=None):
    """Configure axes for plotting a whole genomes's data."""
    pass


def cnv_on_genome(axis, probes, segments, pad, do_trend=False, y_min=None,
                  y_max=None):
    """Plot coverages and CBS calls for all chromosomes on one plot."""
    # Group probes by chromosome (to calculate plotting coordinates)
    if probes:
        chrom_probe_centers = {chrom: 0.5 * (rows['start'] + rows['end'])
                               for chrom, rows in probes.by_chromosome()}
        chrom_sizes = chromosome_sizes(probes)
    else:
        chrom_sizes = chromosome_sizes(segments)

    # Same for segment calls
    chrom_seg_coords = {chrom: list(zip(rows['log2'], rows['start'], rows['end']))
                        for chrom, rows in segments.by_chromosome()
                       } if segments else {}

    x_starts = plot_x_dividers(axis, chrom_sizes, pad)
    x = []
    seg_lines = []  # y-val, x-start, x-end
    for chrom, curr_offset in list(x_starts.items()):
        if probes:
            x.extend(chrom_probe_centers[chrom] + curr_offset)
        if chrom in chrom_seg_coords:
            seg_lines.extend((c[0], c[1] + curr_offset, c[2] + curr_offset)
                             for c in chrom_seg_coords[chrom])

    # Configure axes etc.
    axis.axhline(color='k')
    axis.set_ylabel("Copy ratio (log2)")
    if not (y_min and y_max):
        if segments:
            # Auto-scale y-axis according to segment mean-coverage values
            seg_auto_vals = segments[(segments.chromosome != 'chr6') &
                                     (segments.chromosome != 'chrY')]['log2']
            if not y_min:
                y_min = min(seg_auto_vals.min() - .2, -1.5)
            if not y_max:
                y_max = max(seg_auto_vals.max() + .2, 1.5)
        else:
            if not y_min:
                y_min = -2.5
            if not y_max:
                y_max = 2.5
    axis.set_ylim(y_min, y_max)

    # Plot points
    if probes:
        axis.scatter(x, probes['log2'], color=POINT_COLOR, edgecolor='none',
                     alpha=0.2, marker='.')
        # Add a local trend line
        if do_trend:
            axis.plot(x, _smooth_genome_log2(probes, smoothing.smoothed, 150),
                      color=POINT_COLOR, linewidth=2, zorder=-1)
    # Plot segments
    for seg_line in seg_lines:
        y1, x1, x2 = seg_line
        axis.plot((x1, x2), (y1, y1),
                  color=SEG_COLOR, linewidth=3, solid_capstyle='round')


def _smooth_genome_log2(cnarr, smooth_func, width):
    """Fit a trendline through bin log2 ratios, handling chromosome boundaries.

    Returns
    -------
    np.array
        Smoothed log2 values, calculated with `smooth_func` and `width`, equal
        in length to `cnarr`.
    """
    # ENH: also split by centromeres (long internal gaps -- see PSCBS)
    # ENH: use pandas groupby
    out = {chrom: smooth_func(subcna['log2'], width)
           for chrom, subcna in cnarr.by_chromosome()}
    return np.concatenate(
        [out[chrom] for chrom in sorted(out, key=core.sorter_chrom)])


def snv_on_genome(axis, variants, chrom_sizes, segments, do_trend, pad,
                  do_boost=False):
    """Plot a scatter-plot of SNP chromosomal positions and shifts."""
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("VAF")
    x_starts = plot_x_dividers(axis, chrom_sizes, pad)

    # Calculate the coordinates of plot components
    chrom_snvs = dict(variants.by_chromosome())
    x_posns_chrom = {}
    y_posns_chrom = {}
    trends = []
    for chrom, curr_offset in x_starts.items():
        snvs = chrom_snvs.get(chrom, [])
        if not len(snvs):
            x_posns_chrom[chrom] = []
            y_posns_chrom[chrom] = []
            continue
        posns = np.asfarray(snvs["start"])
        x_posns = posns + curr_offset
        if do_boost:
            vafs = snvs.tumor_boost()
        else:
            vafs = np.asfarray(snvs["alt_freq"])
        x_posns_chrom[chrom] = x_posns
        y_posns_chrom[chrom] = vafs
        # Trend bars: always calculated, only shown on request
        if segments:
            # Draw average VAF within each segment
            for v_start, v_end, v_freq in group_snvs_by_segments(posns, vafs,
                                                                 segments,
                                                                 chrom):
                trends.append((v_start + curr_offset, v_end + curr_offset,
                               v_freq))
        else:
            # Draw chromosome-wide average VAF
            for mask_vaf in ((vafs > 0.5), (vafs <= 0.5)):
                if sum(mask_vaf) > 1:
                    these_posns = x_posns[mask_vaf]
                    trends.append((these_posns[0], these_posns[-1],
                                   np.median(vafs[mask_vaf])))

    # Test for significant shifts in VAF
    # ENH - use segments if provided
    #   if significant, colorize those points / that median line
    sig_chroms = [] # test_loh(partition_by_chrom(chrom_snvs))

    # Render significantly shifted heterozygous regions separately
    x_posns = []
    y_posns = []
    x_posns_sig = []
    y_posns_sig = []
    for chrom in chrom_sizes:
        posns = x_posns_chrom[chrom]
        vafs = y_posns_chrom[chrom]
        if chrom in sig_chroms:
            x_posns_sig.extend(posns)
            y_posns_sig.extend(vafs)
        else:
            x_posns.extend(posns)
            y_posns.extend(vafs)

    # Plot the points
    axis.scatter(x_posns, y_posns, color=POINT_COLOR, edgecolor='none',
                 alpha=0.2, marker='.')
    axis.scatter(x_posns_sig, y_posns_sig, color='salmon', edgecolor='none',
                 alpha=0.3)
    # Add trend lines to each chromosome
    if do_trend or segments:
        # Draw a line across each chromosome at the median shift level
        for x_start, x_end, y_trend in trends:
            # ENH: color by segment gain/loss
            axis.plot([x_start, x_end], [y_trend, y_trend],
                      color='#C0C0C0', linewidth=2, zorder=-1,
                      solid_capstyle='round')



# XXX use by_ranges
def group_snvs_by_segments(snv_posns, snv_freqs, segments, chrom=None):
    """Group SNP allele frequencies by segment.

    Yields
    ------
    tuple
        (start, end, value)
    """
    if chrom:
        segments = segments.select(chromosome=chrom)
    seg_starts = segments.start
    # Assign a segment number to each variant, basically
    indices = np.maximum(seg_starts.searchsorted(snv_posns), 1) - 1
    for i in sorted(set(indices)):
        mask_in_seg = (indices == i)
        freqs = snv_freqs[mask_in_seg]
        posns = snv_posns[mask_in_seg]
        # Separately emit VAFs above and below .5 for plotting
        mask_above_mid = (freqs > 0.5)
        for mask_vaf in (mask_above_mid, ~mask_above_mid):
            if sum(mask_vaf) > 1:
                these_posns = posns[mask_vaf]
                yield (these_posns[0], these_posns[-1],
                       np.median(freqs[mask_vaf]))


def plot_x_dividers(axis, chrom_sizes, pad):
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


def partition_by_chrom(chrom_snvs):
    """Group the tumor shift values by chromosome (for statistical testing)."""
    chromnames = set(chrom_snvs.keys())
    bins = {key: {'thisbin': [], 'otherbins': []}
            for key in chrom_snvs}
    for thischrom, snvs in chrom_snvs.items():
        shiftvals = np.array([abs(v[2]) for v in snvs])
        bins[thischrom]['thisbin'].extend(shiftvals)
        for otherchrom in chromnames:
            if otherchrom == thischrom:
                continue
            bins[otherchrom]['otherbins'].extend(shiftvals)
    return bins


def test_loh(bins, alpha=0.0025):
    """Test each chromosome's SNP shifts and the combined others'.

    The statistical test is Mann-Whitney, a one-sided non-parametric test for
    difference in means.
    """
    # TODO - this doesn't work right if there are many shifted regions
    try:
        from scipy import stats
    except ImportError:
        # SciPy not installed; can't test for significance
        return []

    significant_chroms = []
    for chrom, partitions in bins.items():
        these_shifts = np.array(partitions['thisbin'], np.float_)
        other_shifts = np.array(partitions['otherbins'], np.float_)
        if len(these_shifts) < 20:
            logging.info("Too few points (%d) to test chrom %s",
                         len(these_shifts), chrom)
        elif these_shifts.mean() > other_shifts.mean():
            logging.debug("\nThese ~= %f (N=%d), Other ~= %f (N=%d)",
                          these_shifts.mean(), len(these_shifts),
                          other_shifts.mean(), len(other_shifts))
            u, prob = stats.mannwhitneyu(these_shifts, other_shifts)
            logging.info("Mann-Whitney - %s: u=%s, p=%s", chrom, u, prob)
            if prob < alpha:
                significant_chroms.append(chrom)

    return significant_chroms


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
        orig_names = set()
        for oname in set(gene_probes['gene']):
            orig_names.update(oname.split(','))
        all_coords[chrom][start, end].update(orig_names)
    # Consolidate each region's gene names into a string
    uniq_coords = {}
    for chrom, hits in all_coords.items():
        uniq_coords[chrom] = [(start, end, ",".join(sorted(orig_names)))
                             for (start, end), orig_names in hits.items()]
    return uniq_coords


def gene_coords_by_range(probes, chrom, start, end,
                         ignore=params.IGNORE_GENE_NAMES):
    """Find the chromosomal position of all genes in a range.

    Returns
    -------
    dict
        Of: {chromosome: [(start, end, gene), ...]}
    """
    ignore += ('Background',)
    # Tabulate the genes in the selected region
    genes = collections.OrderedDict()
    for row in probes.in_range(chrom, start, end):
        name = str(row.gene)
        if name in ignore:
            continue
        if name in genes:
            genes[name][1] = row.end
        else:
            genes[name] = [row.start, row.end]
    # Reorganize the data structure
    return {chrom: [(gstart, gend, name)
                    for name, (gstart, gend) in list(genes.items())]}


def unpack_range(a_range):
    """Extract chromosome, start, end from a string or tuple.

    Examples::

        "chr1" -> ("chr1", None, None)
        "chr1:100-123" -> ("chr1", 100, 123)
        ("chr1", 100, 123) -> ("chr1", 100, 123)
    """
    if not a_range:
        return None, None, None
    if isinstance(a_range, basestring):
        if ':' in a_range or '-' in a_range:
            return parse_range_text(a_range)
        return a_range, None, None
    if isinstance(a_range, (list, tuple)) and len(a_range) == 3:
        return tuple(a_range)
    raise ValueError("Not a range: %r" % a_range)


def parse_range_text(text):
    """Parse a chromosomal range specification.

    Parameters
    ----------
    text : string
        Range specification, which should look like ``chr1:1234-5678`` or
        ``chr1:1234-`` or ``chr1:-5678``, where missing start becomes 0 and
        missing end becomes None.
    """
    try:
        chrom, rest = text.split(':')
        start, end = rest.split('-')
        start = int(start) if start else 0
        end = int(end) if end else None
        return chrom, start, end
    except Exception:
        raise ValueError("Invalid range spec: " + text
                         + " (should be like: chr1:2333000-2444000)")
