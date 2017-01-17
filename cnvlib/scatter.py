"""Command-line interface and corresponding API for CNVkit."""
# NB: argparse CLI definitions and API functions are interwoven:
#   "_cmd_*" handles I/O and arguments processing for the command
#   "do_*" runs the command's functionality as an API
from __future__ import absolute_import, division, print_function

import collections
import logging

from matplotlib import pyplot

from . import core, plots
from .cnary import CopyNumArray as CNA


def do_scatter(cnarr, segments=None, variants=None,
               show_range=None, show_gene=None,
               background_marker=None, do_trend=False, window_width=1e6,
               y_min=None, y_max=None, title=None,
               segment_color=plots.SEG_COLOR):
    """Plot probe log2 coverages and segmentation calls together."""
    if not show_gene and not show_range:
        genome_scatter(cnarr, segments, variants, do_trend, y_min, y_max, title,
                       segment_color)
    else:
        chromosome_scatter(cnarr, segments, variants, show_range, show_gene,
                           background_marker, do_trend, window_width, y_min,
                           y_max, title, segment_color)


def genome_scatter(cnarr, segments=None, variants=None, do_trend=False,
                   y_min=None, y_max=None, title=None,
                   segment_color=plots.SEG_COLOR):
    """Plot all chromosomes, concatenated on one plot."""
    if (cnarr or segments) and variants:
        # Lay out top 3/5 for the CN scatter, bottom 2/5 for SNP plot
        axgrid = pyplot.GridSpec(5, 1, hspace=.85)
        axis = pyplot.subplot(axgrid[:3])
        axis2 = pyplot.subplot(axgrid[3:], sharex=axis)
        # Place chromosome labels between the CNR and SNP plots
        axis2.tick_params(labelbottom=False)
        chrom_sizes = plots.chromosome_sizes(cnarr or segments)
        plots.snv_on_genome(axis2, variants, chrom_sizes, segments,
                            do_trend)
    else:
        _fig, axis = pyplot.subplots()
    if title is None:
        title = (cnarr or segments or variants).sample_id
    if cnarr or segments:
        axis.set_title(title)
        plots.cnv_on_genome(axis, cnarr, segments, do_trend, y_min, y_max,
                            segment_color=segment_color)
    else:
        axis.set_title("Variant allele frequencies: %s" % title)
        chrom_sizes = collections.OrderedDict(
            (chrom, subarr["end"].max())
            for chrom, subarr in variants.by_chromosome())
        plots.snv_on_genome(axis, variants, chrom_sizes, segments, do_trend)


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
    chrom, start, end = plots.unpack_range(show_range)
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
            plots.snv_on_chromosome(axis2, sel_snvs, sel_seg, genes, do_trend)
        else:
            _fig, axis = pyplot.subplots()
            axis.set_xlabel("Position (Mb)")
        plots.cnv_on_chromosome(axis, sel_probes, sel_seg, genes,
                                background_marker=background_marker,
                                do_trend=do_trend, x_limits=window_coords,
                                y_min=y_min, y_max=y_max,
                                segment_color=segment_color)
    elif variants:
        # Only plot SNVs in a single-panel layout
        _fig, axis = pyplot.subplots()
        plots.snv_on_chromosome(axis, sel_snvs, sel_seg, genes, do_trend)

    if title is None:
        title = "%s %s" % ((cnarr or segments or variants).sample_id, chrom)
    axis.set_title(title)
