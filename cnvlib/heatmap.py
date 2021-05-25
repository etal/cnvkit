"""The 'heatmap' command."""
import collections
import logging

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from skgenome.rangelabel import unpack_range

from . import plots


def do_heatmap(cnarrs, show_range=None, do_desaturate=False, by_bin=False, vertical=False, ax=None):
    """Plot copy number for multiple samples as a heatmap."""
    if ax is None:
        _fig, axis = plt.subplots()
    else:
        axis = ax

    # List sample names on the y-axis
    if not vertical:
        axis.set_yticks([i + 0.5 for i in range(len(cnarrs))])
        axis.set_yticklabels([c.sample_id for c in cnarrs])
        axis.set_ylim(0, len(cnarrs))
        axis.invert_yaxis()
        axis.set_ylabel("Samples")
    else:
        axis.xaxis.set_major_locator(mticker.FixedLocator([i+1 for i in range(len(cnarrs))]))
        axis.set_xticklabels([c.sample_id for c in cnarrs], rotation=-60)
        axis.set_xlim(0, len(cnarrs))
        #axis.invert_xaxis()
        axis.set_xlabel("Samples")
    
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
        points["log2"] = cna.log2
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
    
    dict_log2 = {}
    if show_range:
        # Lay out only the selected chromosome
        # Set x-axis the chromosomal positions (in Mb), title as the selection
        if by_bin:
            MB = 1
            axis.set_xlabel("Position (bin)")
        else:
            MB = plots.MB
            axis.set_xlabel("Position (Mb)")
        
        if not vertical:
            axis.set_xlim((r_start or 0) * MB,
                          (r_end or chrom_sizes[r_chrom]) * MB)
            print(axis.get_xlim())
            axis.set_title(show_range)
            axis.tick_params(which='both', direction='out')
            axis.get_xaxis().tick_bottom()
            axis.get_yaxis().tick_left()
        else:
            axis.set_ylim((r_start or 0) * MB,
                          (r_end or chrom_sizes[r_chrom]) * MB)
            axis.set_title(show_range)
            axis.tick_params(which='both', direction='out')
            #axis.get_yaxis().tick_bottom() # 'YAxis' object has no attribute 'tick_bottom'
            #axis.get_xaxis().tick_left() # 'XAxis' object has no attribute 'tick_left'
            
        # Plot the individual probe/segment coverages
        for i, sample in enumerate(sample_data):
            sampl_crow = sample[r_chrom]
            if not len(sampl_crow):
                logging.warning("Sample #%d has no datapoints in selection %s",
                                i+1, show_range)
            sampl_crow["start"] *= MB
            sampl_crow["end"] *= MB
            dict_log2[i] = sampl_crow.set_index(['start', 'end']).log2

    else:
        # Lay out chromosome dividers and x-axis labels
        # (Just enough padding to avoid overlap with the divider line)
        chrom_offsets = plots.plot_x_dividers(axis, chrom_sizes, 1, verticaled=vertical)
        # Plot the individual probe/segment coverages
        for i, sample in enumerate(sample_data):
            all_crows = []
            for chrom, curr_offset in chrom_offsets.items():
                crow = sample[chrom]
                if len(crow):
                    crow["start"] += curr_offset
                    crow["end"] += curr_offset
                else:
                    logging.warning("Sample #%d has no datapoints", i+1)
                all_crows.append(crow)
            
            sampl_crow = pd.concat(all_crows, axis='index')
            dict_log2[i] = sampl_crow.set_index(['start', 'end']).log2

    log2_df = pd.DataFrame.from_dict(dict_log2)
    # Need to explicitly insert NaN-filled rows in-between 2 discontiguous intervals
    log2_df.reset_index(inplace=True)
    compt = 0
    for i in range(1, len(log2_df.index)):
        end_previous = log2_df.iloc[i-1].end ; start_current = log2_df.iloc[i].start
        if end_previous != start_current: # Discontinous
            compt += 1
            log2_df.loc[i-1+0.5, :] = [end_previous, start_current] + [np.nan] * len(cnarrs)
    log2_df = log2_df.sort_index().set_index(['start', 'end'])
    print("INSERTED", compt, "EMTPY intervals (log2=NaN for all samples)")
    
    start_val = list(log2_df.index.get_level_values('start'))
    end_val = list(log2_df.index.get_level_values('end'))
    # "If shading='flat' (which is default) the dimensions of X and Y should be one greater than those of C":
    start2plt = np.array(start_val + [end_val[-1]])
    sampl2plt = np.array(range(len(cnarrs) + 1))
    
    if not vertical: # BEWARE 'normal old view' == 'pcolor on transposed_df' 
        fixed_start = start_val
        fixed_end = end_val
        dat2plot = log2_df.transpose()
        X_pcolor, Y_pcolor = start2plt, sampl2plt
    else:
        fixed_start = np.array(log2_df.columns)
        fixed_end = np.array(log2_df.columns + 1)
        dat2plot = log2_df
        X_pcolor, Y_pcolor = sampl2plt, start2plt # INVERSION COMPARED TO 'not_vertical'

    cMap = plt.get_cmap('bwr')
    # 'CenteredNorm' looks like 'desaturate' process
    # if do_desaturate and hasattr(mpl.colors, 'CenteredNorm'): # Requires matplotlib >= 3.4.2...
    if False: # NO correct norm yet for 'desaturate'
        norm = mpl.colors.CenteredNorm(halfrange=5) # 'halfrange=5' is empirically a good value
        im = axis.pcolormesh(X_pcolor, Y_pcolor, dat2plot, norm=norm, cmap=cMap)
    else:
        im = axis.pcolormesh(X_pcolor, Y_pcolor, dat2plot, vmin=-1.33, vmax=1.33, cmap=cMap)
    cbar = plt.colorbar(im, ax=axis, fraction=0.04, pad=0.03, shrink=0.6)
    cbar.set_label("log2", labelpad=0)

    return axis

