"""The 'heatmap' command."""
from __future__ import absolute_import, division, print_function
from builtins import zip

import collections
import logging

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from skgenome.rangelabel import unpack_range

from . import plots


def do_heatmap(cnarrs, show_range=None, do_desaturate=False, by_bin=False):
    """Plot copy number for multiple samples as a heatmap."""
    _fig, axis = plt.subplots()
    set_colorbar(axis)

    # List sample names on the y-axis
    axis.set_yticks([i + 0.5 for i in range(len(cnarrs))])
    axis.set_yticklabels([c.sample_id for c in cnarrs])
    axis.set_ylim(0, len(cnarrs))
    axis.invert_yaxis()
    axis.set_ylabel("Samples")
    if hasattr(axis, 'set_facecolor'):
        # matplotlib >= 2.0
        axis.set_facecolor('#DDDDDD')
    else:
        # Older matplotlib
        axis.set_axis_bgcolor('#DDDDDD')

    if by_bin and show_range:
        try:
            a_cnarr = next(c for c in cnarrs if "probes" not in c)
        except StopIteration:
            r_chrom, r_start, r_end = unpack_range(show_range)
            if r_start is not None or r_end is not None:
                raise ValueError("Need at least 1 .cnr input file if --by-bin "
                                 "(by_bin) and --chromosome (show_range) are "
                                 "both used to specify a sub-chromosomal "
                                 "region.")
        else:
            logging.info("Using sample %s to map %s to bin coordinates",
                         a_cnarr.sample_id, show_range)
            r_chrom, r_start, r_end = plots.translate_region_to_bins(show_range,
                                                                     a_cnarr)
    else:
        r_chrom, r_start, r_end = unpack_range(show_range)
    if r_start is not None or r_end is not None:
        logging.info("Showing log2 ratios in range %s:%d-%s",
                     r_chrom, r_start or 0, r_end or '*')
    elif r_chrom:
        logging.info("Showing log2 ratios on chromosome %s", r_chrom)

    # Closes over do_desaturate
    def cna2df(cna):
        """Extract a dataframe of plotting points from a CopyNumArray."""
        points = cna.data.loc[:, ["start", "end"]]
        points["color"] = cna.log2.apply(plots.cvg2rgb, args=(do_desaturate,))
        return points

    # Group each file's probes/segments by chromosome
    sample_data = [collections.defaultdict(list) for _c in cnarrs]
    # Calculate the size (max endpoint value) of each chromosome
    chrom_sizes = collections.OrderedDict()
    for i, cnarr in enumerate(cnarrs):
        if by_bin:
            cnarr = plots.update_binwise_positions_simple(cnarr)

        if r_chrom:
            subcna = cnarr.in_range(r_chrom, r_start, r_end, mode="trim")
            sample_data[i][r_chrom] = cna2df(subcna)
            chrom_sizes[r_chrom] = max(subcna.end.iat[-1] if subcna else 0,
                                       chrom_sizes.get(r_chrom, 0))
        else:
            for chrom, subcna in cnarr.by_chromosome():
                sample_data[i][chrom] = cna2df(subcna)
                chrom_sizes[chrom] = max(subcna.end.iat[-1] if subcna else 0,
                                         chrom_sizes.get(r_chrom, 0))

    # Closes over axis
    def plot_sample_chrom(i, sample):
        """Draw the given coordinates and colors as a horizontal series."""
        xranges = [(start, end - start)
                   for start, end in zip(sample.start, sample.end)]
        bars = BrokenBarHCollection(xranges, (i, i+1),
                                    edgecolors="none",
                                    facecolors=sample["color"])
        axis.add_collection(bars)

    if show_range:
        # Lay out only the selected chromosome
        # Set x-axis the chromosomal positions (in Mb), title as the selection
        if by_bin:
            MB = 1
            axis.set_xlabel("Position (bin)")
        else:
            MB = plots.MB
            axis.set_xlabel("Position (Mb)")
        axis.set_xlim((r_start or 0) * MB,
                      (r_end or chrom_sizes[r_chrom]) * MB)
        axis.set_title(show_range)
        axis.tick_params(which='both', direction='out')
        axis.get_xaxis().tick_bottom()
        axis.get_yaxis().tick_left()
        # Plot the individual probe/segment coverages
        for i, sample in enumerate(sample_data):
            crow = sample[r_chrom]
            if not len(crow):
                logging.warning("Sample #%d has no datapoints in selection %s",
                                i+1, show_range)
            crow["start"] *= MB
            crow["end"] *= MB
            plot_sample_chrom(i, crow)

    else:
        # Lay out chromosome dividers and x-axis labels
        # (Just enough padding to avoid overlap with the divider line)
        chrom_offsets = plots.plot_x_dividers(axis, chrom_sizes, 1)
        # Plot the individual probe/segment coverages
        for i, sample in enumerate(sample_data):
            for chrom, curr_offset in chrom_offsets.items():
                crow = sample[chrom]
                if len(crow):
                    crow["start"] += curr_offset
                    crow["end"] += curr_offset
                    plot_sample_chrom(i, crow)
                else:
                    logging.warning("Sample #%d has no datapoints", i+1)

    return axis


def set_colorbar(axis):
    # Create our colormap
    # ENH: refactor to use colormap to colorize the BrokenBarHCollection
    #   - maybe also refactor plots.cvg2rgb
    cmap = mpl.colors.LinearSegmentedColormap.from_list('cnvheat',
        [(0, 0, .75),
         (1, 1, 1),
         (.75, 0, 0)])
    # Add a colorbar
    norm = mpl.colors.Normalize(-1.33, 1.33)
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(np.linspace(-1.33, 1.33, 30))
    cbar = plt.colorbar(mappable, ax=axis, orientation='vertical',
                        fraction=0.04, pad=0.03, shrink=0.6,
                        #  label="log2",
                        ticks=(-1, 0, 1))
    cbar.set_label("log2", labelpad=0)
