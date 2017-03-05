"""Command-line interface and corresponding API for CNVkit."""
# NB: argparse CLI definitions and API functions are interwoven:
#   "_cmd_*" handles I/O and arguments processing for the command
#   "do_*" runs the command's functionality as an API
from __future__ import absolute_import, division, print_function

import collections
import logging

import numpy as np
from matplotlib import pyplot
from skgenome.rangelabel import unpack_range

from . import core, plots, smoothing
from .plots import MB
from .cnary import CopyNumArray as CNA

HIGHLIGHT_COLOR = 'gold'
POINT_COLOR = '#606060'
SEG_COLOR = 'darkorange'
TREND_COLOR = '#C0C0C0'


def do_scatter(cnarr, segments=None, variants=None,
               show_range=None, show_gene=None,
               background_marker=None, do_trend=False, window_width=1e6,
               y_min=None, y_max=None, title=None, segment_color=SEG_COLOR):
    """Plot probe log2 coverages and segmentation calls together."""
    if not show_gene and not show_range:
        genome_scatter(cnarr, segments, variants, do_trend, y_min, y_max, title,
                       segment_color)
    else:
        chromosome_scatter(cnarr, segments, variants, show_range, show_gene,
                           background_marker, do_trend, window_width, y_min,
                           y_max, title, segment_color)


# === Genome-level scatter plots ===

def genome_scatter(cnarr, segments=None, variants=None, do_trend=False,
                   y_min=None, y_max=None, title=None, segment_color=SEG_COLOR):
    """Plot all chromosomes, concatenated on one plot."""
    if (cnarr or segments) and variants:
        # Lay out top 3/5 for the CN scatter, bottom 2/5 for SNP plot
        axgrid = pyplot.GridSpec(5, 1, hspace=.85)
        axis = pyplot.subplot(axgrid[:3])
        axis2 = pyplot.subplot(axgrid[3:], sharex=axis)
        # Place chromosome labels between the CNR and SNP plots
        axis2.tick_params(labelbottom=False)
        chrom_sizes = plots.chromosome_sizes(cnarr or segments)
        snv_on_genome(axis2, variants, chrom_sizes, segments, do_trend)
    else:
        _fig, axis = pyplot.subplots()
    if title is None:
        title = (cnarr or segments or variants).sample_id
    if cnarr or segments:
        axis.set_title(title)
        cnv_on_genome(axis, cnarr, segments, do_trend, y_min, y_max,
                      segment_color=segment_color)
    else:
        axis.set_title("Variant allele frequencies: %s" % title)
        chrom_sizes = collections.OrderedDict(
            (chrom, subarr["end"].max())
            for chrom, subarr in variants.by_chromosome())
        snv_on_genome(axis, variants, chrom_sizes, segments, do_trend)



def cnv_on_genome(axis, probes, segments, do_trend=False, y_min=None,
                  y_max=None, segment_color=SEG_COLOR):
    """Plot coverages and CBS calls for all chromosomes on one plot."""
    # Group probes by chromosome (to calculate plotting coordinates)
    if probes:
        chrom_probe_centers = {chrom: 0.5 * (rows['start'] + rows['end'])
                               for chrom, rows in probes.by_chromosome()}
        chrom_sizes = plots.chromosome_sizes(probes)
    else:
        chrom_sizes = plots.chromosome_sizes(segments)

    # Same for segment calls
    chrom_seg_coords = {chrom: list(zip(rows['log2'], rows['start'], rows['end']))
                        for chrom, rows in segments.by_chromosome()
                       } if segments else {}

    x_starts = plots.plot_x_dividers(axis, chrom_sizes)
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
                  color=segment_color, linewidth=3, solid_capstyle='round')


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
    out = [smooth_func(subcna['log2'], width)
           for _chrom, subcna in cnarr.by_chromosome()]
    return np.concatenate(out)


def snv_on_genome(axis, variants, chrom_sizes, segments, do_trend):
    """Plot a scatter-plot of SNP chromosomal positions and shifts."""
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("VAF")
    x_starts = plots.plot_x_dividers(axis, chrom_sizes)

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
        posns = snvs['start'].values
        x_posns = posns + curr_offset
        vafs = snvs['alt_freq'].values
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
                      color=TREND_COLOR, linewidth=2, zorder=-1,
                      solid_capstyle='round')


# === Chromosome-level scatter plots ===

def chromosome_scatter(cnarr, segments, variants, show_range, show_gene,
                       background_marker, do_trend, window_width, y_min, y_max,
                       title, segment_color):
    """Plot a specified region on one chromosome.

    Possibilities::

             Options | Shown
        ------------ | --------
        -c      | -g | Genes | Region
        ------- | -- | ----- | ------
        -       | +  | given | auto: gene(s) + margin
        chr     | -  | none  | whole chrom
        chr     | +  | given | whole chrom
        chr:s-e | -  | all   | given
        chr:s-e | +  | given | given

    """
    chrom, start, end = unpack_range(show_range)
    window_coords = ()
    genes = []
    if show_gene:
        gene_names = show_gene.split(',')
        # Scan for probes matching the specified gene
        gene_coords = plots.gene_coords_by_name(cnarr or segments,
                                                gene_names)
        if len(gene_coords) != 1:
            raise ValueError("Genes %s are split across chromosomes %s"
                             % (show_gene, list(gene_coords.keys())))
        g_chrom, genes = gene_coords.popitem()
        if chrom:
            # Confirm that the selected chromosomes match
            core.assert_equal("Chromosome also selected by region (-c) "
                              "does not match",
                              **{"chromosome": chrom,
                                 "gene(s)": g_chrom})
        else:
            chrom = g_chrom
        # Set the display window to the selected genes +/- a margin
        genes.sort()
        window_coords = (max(0, genes[0][0] - window_width),
                         genes[-1][1] + window_width)

    if start is not None or end is not None:
        # Default selection endpoint to the maximum chromosome position
        if not end:
            end = (cnarr or segments or variants
                  ).filter(chromosome=chrom).end.iat[-1]
        if window_coords:
            # Genes were specified, & window was set around them
            if start > window_coords[0] or end < window_coords[1]:
                raise ValueError("Selected gene region " + chrom +
                                 (":%d-%d" % window_coords) +
                                 " is outside specified region " +
                                 show_range)
        window_coords = (max(0, start - window_width), end + window_width)
        if cnarr and not genes:
            genes = plots.gene_coords_by_range(cnarr, chrom, start, end)[chrom]
        if not genes and window_width > (end - start) / 10.0:
            # No genes in the selected region, so highlight the region
            # itself (unless the selection is ~entire displayed window)
            logging.info("No genes found in selection; will show the "
                         "selected region itself instead")
            genes = [(start, end, "Selection")]
    elif show_range and window_coords:
        # Specified range is only chrom, no start-end
        # Reset window around selected genes to show the whole chromosome
        window_coords = ()

    # Prune plotted elements to the selected region
    sel_probes = (cnarr.in_range(chrom, *window_coords)
                  if cnarr else CNA([]))
    sel_seg = (segments.in_range(chrom, *window_coords, mode='trim')
               if segments else CNA([]))
    sel_snvs = (variants.in_range(chrom, *window_coords)
                if variants else None)
    logging.info("Showing %d probes and %d selected genes in region %s",
                 len(sel_probes), len(genes),
                 (chrom + ":%d-%d" % window_coords if window_coords else chrom))

    # Create plots
    if cnarr or segments:
        # Plot CNVs at chromosome level
        if variants:
            # Lay out top 3/5 for the CN scatter, bottom 2/5 for SNP plot
            axgrid = pyplot.GridSpec(5, 1, hspace=.5)
            axis = pyplot.subplot(axgrid[:3])
            axis2 = pyplot.subplot(axgrid[3:], sharex=axis)
            # Plot allele freqs for only the selected region
            snv_on_chromosome(axis2, sel_snvs, sel_seg, genes, do_trend)
        else:
            _fig, axis = pyplot.subplots()
            axis.set_xlabel("Position (Mb)")
        cnv_on_chromosome(axis, sel_probes, sel_seg, genes,
                          background_marker=background_marker,
                          do_trend=do_trend, x_limits=window_coords,
                          y_min=y_min, y_max=y_max, segment_color=segment_color)
    elif variants:
        # Only plot SNVs in a single-panel layout
        _fig, axis = pyplot.subplots()
        snv_on_chromosome(axis, sel_snvs, sel_seg, genes, do_trend)

    if title is None:
        title = "%s %s" % ((cnarr or segments or variants).sample_id, chrom)
    axis.set_title(title)


def cnv_on_chromosome(axis, probes, segments, genes, background_marker=None,
                      do_trend=False, x_limits=None, y_min=None, y_max=None,
                      segment_color=SEG_COLOR):
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
    if x_limits:
        x_min, x_max = x_limits
        axis.set_xlim(x_min * MB, x_max * MB)
    else:
        set_xlim_from(axis, probes, segments)
    setup_chromosome(axis, y_min, y_max, "Copy ratio (log2)")
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
                      color=segment_color, linewidth=4, solid_capstyle='round')


def snv_on_chromosome(axis, variants, segments, genes, do_trend):
    # TODO set x-limits if not already done for probes/segments
    # set_xlim_from(axis, None, segments, variants)
    # setup_chromosome(axis, 0.0, 1.0, "VAF")
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("VAF")
    axis.set_xlabel("Position (Mb)")
    axis.get_yaxis().tick_left()
    axis.get_xaxis().tick_top()
    axis.tick_params(which='both', direction='out',
                     labelbottom=False, labeltop=False)

    x_mb = variants['start'] * MB
    y = variants['alt_freq'].values
    axis.scatter(x_mb, y, color=POINT_COLOR, alpha=0.3)
    # TODO - highlight genes/selection
    if segments:
        # Draw average VAF within each segment
        # TODO coordinate this w/ do_trend
        posns = variants['start'].values # * MB
        for v_start, v_end, v_freq in group_snvs_by_segments(posns, y,
                                                             segments):
            # ENH: color by segment gain/loss
            axis.plot([v_start * MB, v_end * MB], [v_freq, v_freq],
                      color='#C0C0C0', linewidth=2, #zorder=1,
                      solid_capstyle='round')


def set_xlim_from(axis, probes=None, segments=None, variants=None):
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
    if max_x <= min_x:
        if min_x != np.inf:
            logging.warn("*WARNING* selection start %s > end %s; did you "
                         "correctly sort the input file by genomic location?",
                         min_x, max_x)
        raise ValueError("No usable data points to plot out of "
                         "%d probes, %d segments, %d variants"
                         % (len(probes) if probes else 0,
                            len(segments) if segments else 0,
                            len(variants) if variants else 0))
    axis.set_xlim(min_x * MB, max_x * MB)


def setup_chromosome(axis, y_min=None, y_max=None, y_label=None):
    """Configure axes for plotting a single chromosome's data."""
    if y_min and y_max:
        axis.set_ylim(y_min, y_max)
        if y_min < 0 < y_max:
            axis.axhline(color='k')
    if y_label:
        axis.set_ylabel(y_label)
    axis.tick_params(which='both', direction='out')
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()


# === Shared ===

# XXX use by_ranges
def group_snvs_by_segments(snv_posns, snv_freqs, segments, chrom=None):
    """Group SNP allele frequencies by segment.

    Yields
    ------
    tuple
        (start, end, value)
    """
    if chrom:
        segments = segments.filter(chromosome=chrom)
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
