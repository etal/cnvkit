"""Plotting utilities."""
from __future__ import absolute_import, division
from builtins import str
from past.builtins import basestring

import collections
import logging
import math

import numpy as np

from . import core, params


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
