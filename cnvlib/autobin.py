"""Estimate reasonable bin sizes from BAM read counts or depths."""
from __future__ import absolute_import, division, print_function

import tempfile

import numpy as np
import pandas as pd
import pysam

from . import coverage, tabio
from .descriptives import weighted_median


def average_depth(rc_table, read_length):
    """Estimate the average read depth across the genome.

    Returns
    -------
    float
        Median of the per-chromosome mean read depths, weighted by chromosome
        size.
    """
    mean_depths = read_length * rc_table.mapped / rc_table.length
    return weighted_median(mean_depths, rc_table.length)


def update_chrom_length(rc_table, regions):
    if regions is not None and len(regions):
        chrom_sizes = region_size_by_chrom(regions)
        rc_table = rc_table.merge(chrom_sizes, on='chromosome', how='inner')
        rc_table['length'] = rc_table['length_y'] # ?
        rc_table = rc_table.drop(['length_x', 'length_y'], axis=1)
    return rc_table


def region_size_by_chrom(regions):
    chromgroups = regions.data.groupby('chromosome', sort=False)
    # sizes = chromgroups.apply(total_region_size) # XXX
    sizes = [total_region_size(g) for _key, g in chromgroups]
    return pd.DataFrame({'chromosome': regions.chromosome.drop_duplicates(),
                         'length': sizes})


def total_region_size(regions):
    """Aggregate area of all genomic ranges in `regions`."""
    return (regions.end - regions.start).sum()


def sample_region_cov(bam_fname, regions, max_num=100):
    """Calculate read depth in a randomly sampled subset of regions."""
    midsize_regions = sample_midsize_regions(regions, max_num)
    with tempfile.NamedTemporaryFile(suffix='.bed', # mode='w+t'
                                    ) as f:
        tabio.write(regions.as_dataframe(midsize_regions), f, 'bed4')
        f.flush()
        table = coverage.bedcov(f.name, bam_fname, 0)
    table = table[(table.end > table.start) & (table.basecount > 0)]  # Safety
    depths = table.basecount / (table.end - table.start)
    return depths.median()


def sample_midsize_regions(regions, max_num):
    """Randomly select a subset of up to `max_num` regions."""
    sizes = regions.end - regions.start
    lo_size, hi_size = np.percentile(sizes, [25, 75])
    midsize_regions = regions.data[(sizes >= lo_size) & (sizes <= hi_size)]
    if len(midsize_regions) > max_num:
        midsize_regions = midsize_regions.sample(max_num)
    return midsize_regions
