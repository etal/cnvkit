"""Plotting utilities."""
from __future__ import absolute_import, division

import collections
import math
import sys

import numpy
# from matplotlib import pyplot
# pyplot.ioff()

from Bio._py3k import map, range, zip
iteritems = (dict.iteritems if sys.version < 3
             else dict.items)

from . import core, smoothing
from .ngfrills import echo


MB = 1e-6


def plot_genome(axis, probes, segments, pad, do_trend=False):
    """Plot coverages and CBS calls for all chromosomes on one plot."""
    # Group probes by chromosome (to calculate plotting coordinates)
    chrom_probe_centers = {chrom: 0.5 * (rows['start'] + rows['end'])
                           for chrom, rows in probes.by_chromosome()}
    # Same for segment calls
    chrom_seg_coords = {chrom: zip(rows['coverage'], rows['start'], rows['end'])
                        for chrom, rows in segments.by_chromosome()
                       } if segments else []

    chrom_sizes = chromosome_sizes(probes)
    x_starts = plot_x_dividers(axis, chrom_sizes, pad)
    x = []
    seg_lines = []  # y-val, x-start, x-end
    for chrom, curr_offset in x_starts.items():
        x.extend(chrom_probe_centers[chrom] + curr_offset)
        if chrom in chrom_seg_coords:
            seg_lines.extend((c[0], c[1] + curr_offset, c[2] + curr_offset)
                             for c in chrom_seg_coords[chrom])

    # Configure axes etc.
    axis.axhline(color='k')
    axis.set_ylabel("Copy ratio (log2)")
    if segments:
        # Auto-scale y-axis according to segment mean-coverage values
        seg_auto_vals = segments[(segments.chromosome != 'chr6') &
                                 (segments.chromosome != 'chrY')]['coverage']
        axis.set_ylim(min(seg_auto_vals.min() - .2, -1.5),
                      max(seg_auto_vals.max() + .2, 1.5))
    else:
        axis.set_ylim(-2.5, 2.5)

    # Plot points
    axis.scatter(x, probes.coverage, color='#808080', edgecolor='none',
                 alpha=0.15, marker='.')
    # Add a local trend line
    if do_trend:
        axis.plot(x, smoothing.smooth_genome_coverages(probes,
                                                       smoothing.smoothed,
                                                       250),
                  color='#808080', linewidth=1, zorder=-1)
    # Plot segments
    for seg_line in seg_lines:
        axis.hlines(*seg_line, colors='red', linewidth=2)


def plot_chromosome(axis, probes, segments, chromosome, sample, genes,
                    background_marker=None, do_trend=False):
    """Draw a scatter plot of probe values with CBS calls overlaid.

    Argument 'genes' is a list of tuples: (start, end, gene name)
    """
    # Get scatter plot coordinates
    sel_probes = probes[probes['chromosome'] == chromosome]
    x = [probe_center(row) * MB for row in sel_probes]
    y = sel_probes['coverage']
    if 'weight' in sel_probes.dtype.fields:
        w = 46 * sel_probes['weight'] ** 2 + 2
    else:
        w = numpy.repeat(30, len(x))
    is_bg = (sel_probes['gene'] == 'Background')

    # Configure axes
    axis.axhline(color='k')
    axis.set_xlim(max(0, min(x)), max(x))
    axis.set_ylim(limit(min(y) - .1, -5.0, -.3),
                  limit(max(y) + (.25 if genes else .1), .3, 5.0))
    axis.set_ylabel("Copy ratio (log2)")
    axis.set_title("%s %s" % (sample, chromosome))
    axis.tick_params(which='both', direction='out')
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
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
            # Rescale positions from bases to megabases
            gene_start *= MB
            gene_end *= MB
            # Highlight and label gene region
            axis.axvspan(gene_start, gene_end, alpha=0.5, color='gold',
                         zorder=-1)
            axis.text(0.5 * (gene_start + gene_end), min(max(y) + .1, 2.4),
                      gene_name, horizontalalignment='center',
                      rotation=text_rot,
                      size=text_size)
                      # size='small')

    if background_marker in (None, 'o'):
        # Plot targets and antitargets with the same marker
        axis.scatter(x, y, w, color='#808080', alpha=0.4, marker='o')
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
        axis.scatter(x_fg, y_fg, w_fg, color='#808080', alpha=0.4, marker='o')
        axis.scatter(x_bg, y_bg, color='#808080', alpha=0.5,
                     marker=background_marker)

    # Add a local trend line
    if do_trend:
        axis.plot(x, smoothing.smoothed(y, 100),
                    color='#808080', linewidth=1, zorder=-1)

    # Get coordinates for CBS lines & draw them
    if segments:
        for row in segments[segments['chromosome'] == chromosome]:
            axis.hlines(row['coverage'], row['start'] * MB, row['end'] * MB,
                        colors='red', linewidth=2)


def plot_loh(axis, chrom_snvs, chrom_sizes, do_trend, pad):
    """Plot a scatter-plot of SNP chromosomal positions and shifts."""
    axis.set_ylim(0.5, 1.0)
    axis.set_ylabel("VAF")
    x_starts = plot_x_dividers(axis, chrom_sizes, pad)

    # Calculate the coordinates of plot components
    x_posns_chrom = {}
    y_posns_chrom = {}
    trends = []
    for chrom, curr_offset in iteritems(x_starts):
        snvs = chrom_snvs[chrom]
        if not len(snvs):
            x_posns_chrom[chrom] = []
            y_posns_chrom[chrom] = []
            continue
        posns = numpy.array([v[0] for v in snvs], numpy.float_)
        x_posns = posns + curr_offset
        vafs = numpy.array([abs(v[2] - .5) + 0.5 for v in snvs], numpy.float_)
        x_posns_chrom[chrom] = x_posns
        y_posns_chrom[chrom] = vafs
        # Trend bars: always calculated, only shown on request
        trends.append((x_posns[0], x_posns[-1], numpy.median(vafs)))

    # Test for significant shifts in VAF
    # ENH - use segments if provided
    sig_chroms = test_loh(partition_by_chrom(chrom_snvs))
    if sig_chroms:
        echo("Significant LOH shift in chromosomes", ", ".join(sig_chroms))
    else:
        echo("No significant LOH shift in any chromosomes")

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
    axis.scatter(x_posns, y_posns, color='#808080', edgecolor='none', alpha=0.2)
    axis.scatter(x_posns_sig, y_posns_sig, color='salmon', edgecolor='none',
                 alpha=0.3)
    # Add trend lines to each chromosome
    if do_trend:
        # Draw a line across each chromosome at the median shift level
        for x_start, x_end, y_trend in trends:
            axis.plot([x_start, x_end], [y_trend, y_trend],
                      color='#C0C0C0', linewidth=3, zorder=-1)


def plot_x_dividers(axis, chromosome_sizes, pad):
    """Plot vertical dividers and x-axis labels given the chromosome sizes.

    Returns a table of the x-position offsets of each chromosome.

    Draws vertical black lines between each chromosome, with padding.
    Labels each chromosome range with the chromosome name, centered in the
    region, under a tick.
    Sets the x-axis limits to the covered range.
    """
    assert isinstance(chromosome_sizes, collections.OrderedDict)

    x_dividers = []
    x_centers = []
    x_starts = collections.OrderedDict()
    curr_offset = pad
    for label, size in chromosome_sizes.items():
        x_starts[label] = curr_offset
        x_centers.append(curr_offset + 0.5 * size)
        x_dividers.append(curr_offset + size + pad)
        curr_offset += size + 2 * pad

    axis.set_xlim(0, curr_offset)
    for xposn in x_dividers[:-1]:
        axis.axvline(x=xposn, color='k')
    # Use chromosome names as x-axis labels (instead of base positions)
    axis.set_xticks(x_centers)
    axis.set_xticklabels(chromosome_sizes.keys(), rotation=60)
    axis.tick_params(labelsize='small')
    axis.tick_params(axis='x', length=0)
    # axis.tick_params(direction='out') # XXX Y U CRASH ?!
    # axis.get_xaxis().tick_bottom()
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



def limit(x, lower, upper):
    """Limit x to between lower and upper bounds."""
    assert lower < upper
    return max(lower, min(x, upper))


def probe_center(row):
    """Return the midpoint of the probe location."""
    return 0.5 * (row['start'] + row['end'])


def partition_by_chrom(chrom_snvs):
    """Group the tumor shift values by chromosome (for statistical testing)."""
    chromnames = set(chrom_snvs.keys())
    bins = {key: {'thisbin': [], 'otherbins': []}
            for key in chrom_snvs}
    for thischrom, snvs in iteritems(chrom_snvs):
        shiftvals = numpy.array([abs(v[2]) for v in snvs])
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
    # maybe: sort bins by mean, and work downward/upward in order
    # or: start with kruskal-wallis
    try:
        from scipy import stats
    except ImportError:
        # SciPy not installed; can't test for significance
        return []

    significant_chroms = []
    for chrom, partitions in iteritems(bins):
        these_shifts = numpy.array(partitions['thisbin'], numpy.float_)
        other_shifts = numpy.array(partitions['otherbins'], numpy.float_)
        if len(these_shifts) < 20:
            echo("Too few points (%d) to test chrom %s"
                 % (len(these_shifts), chrom))
        elif these_shifts.mean() > other_shifts.mean():
            # DBG
            echo("\nThese ~= %f (N=%d), Other ~= %f (N=%d)" %
                   (these_shifts.mean(), len(these_shifts),
                    other_shifts.mean(), len(other_shifts)))
            # ---
            u, prob = stats.mannwhitneyu(these_shifts, other_shifts)
            echo("Mann-Whitney - %s: u=%s, p=%s" % (chrom, u, prob))
            if prob < alpha:
                significant_chroms.append(chrom)

    return significant_chroms


# ________________________________________
# Utilies used by other modules


def cvg2rgb(cvg, desaturate):
    """Choose a shade of red or blue representing log2-coverage value."""
    cutoff = 1.25  # Coverages above this point are shown with max intensity
    x = min(abs(cvg) / cutoff, 1.0)
    if desaturate:
        # Adjust intensity sigmoidally -- reduce near 0, boost near 1
        # Exponent <1 shifts the fixed point leftward (from x=0.5)
        x = ((1. - math.cos(x * math.pi)) / 2.) ** 0.9
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
def gene_coords_by_name(probes, names):
    """Find the chromosomal position of each named gene in probes.

    Returns a dict: {chromosome: [(start, end, gene name), ...]}
    """
    coords = collections.defaultdict(list)
    for name in names:
        gene_probes = probes[probes.gene == name]
        if not len(gene_probes):
            raise ValueError("No targeted gene named '%s' found" % name)
        # Find the genomic range of this gene's probes
        start = gene_probes['start'].min()
        end = gene_probes['end'].max()
        chrom = core.check_unique(gene_probes['chromosome'], name)
        coords[chrom].append((start, end, name))
    return coords


# XXX should this be a CopyNumArray method?
def gene_coords_by_range(probes, chrom, start, end,
                         skip=('Background', 'CGH', '-')):
    """Find the chromosomal position of all genes in a range.

    Returns a dict: {chromosome: [(start, end, gene), ...]}
    """
    # Tabulate the genes in the selected region
    genes = collections.OrderedDict()
    for row in probes.in_range(chrom, start, end):
        name = str(row['gene'])
        if name in skip:
            continue
        if name in genes:
            genes[name][1] = row['end']
        else:
            genes[name] = [row['start'], row['end']]
    # Reorganize the data structure
    return {chrom: [(start, end, name)
                    for name, (start, end) in genes.items()]}


# XXX not really specific to plots...
def parse_range(text):
    """Parse a chromosomal range specification.

    Range spec string should look like: 'chr1:1234-5678'
    """
    try:
        chrom, rest = text.split(':')
        start, end = map(int, rest.split('-'))
        return chrom, start, end
    except Exception:
        raise ValueError("Invalid range spec: " + text
                         + " (should be like: chr1:2333000-2444000)")

