"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division, print_function
from builtins import map

import logging
import re

from . import tabio
from .params import INSERT_SIZE
from .genome import GenomicArray as GA


def get_background(target_bed, access_bed, avg_bin_size, min_bin_size):
    """Generate background intervals from target intervals.

    Procedure:

    - Invert target intervals
    - Subtract the inverted targets from accessible regions
    - For each of the resulting regions:

        - Shrink by a fixed margin on each end
        - If it's smaller than min_bin_size, skip
        - Divide into equal-size (region_size/avg_bin_size) portions
        - Emit the (chrom, start, end) coords of each portion
    """
    targets = tabio.read_auto(target_bed)
    target_chroms = set(targets.chromosome.unique())
    if access_bed:
        # Chromosomes' accessible sequence regions are given -- use them
        accessible = tabio.read_auto(access_bed)
        access_chroms = set(accessible.chromosome.unique())
        if access_chroms and access_chroms.isdisjoint(target_chroms):
            raise ValueError("Chromosome names in the accessible regions file "
                             "%s %r do not match those in targets %s %r"
                             % (access_bed, tuple(sorted(access_chroms)[:3]),
                                target_bed, tuple(sorted(target_chroms)[:3])))
        # But filter out untargeted alternative contigs and mitochondria
        untgt_chroms = access_chroms - target_chroms
        # Autosomes typically have numeric names, allosomes are X and Y
        is_canonical = re.compile(r"(chr)?(\d+|[XYxy])$")
        if any(is_canonical.match(c) for c in target_chroms):
            chroms_to_skip = [c for c in untgt_chroms
                              if not is_canonical.match(c)]
        else:
            # Alternative contigs have longer names -- skip them
            max_tgt_chr_name_len = max(map(len, target_chroms))
            chroms_to_skip = [c for c in untgt_chroms
                              if len(c) > max_tgt_chr_name_len]
        if chroms_to_skip:
            logging.info("Skipping untargeted chromosomes %s",
                         ' '.join(sorted(chroms_to_skip)))
            skip_idx = accessible.chromosome.isin(chroms_to_skip)
            accessible = accessible[~skip_idx]
    else:
        # Chromosome accessible sequence regions not known -- use heuristics
        # (chromosome length is endpoint of last probe; skip initial
        # <magic number> of bases that are probably telomeric)
        TELOMERE_SIZE = 150000
        accessible = guess_chromosome_regions(targets, TELOMERE_SIZE)

    pad_size = 2 * INSERT_SIZE
    bg_arr = (accessible.resize_ranges(-pad_size)
              .subtract(targets.resize_ranges(pad_size)))
    bg_arr = bg_arr.subdivide(avg_bin_size, min_bin_size)
    bg_arr['gene'] = 'Background'
    return bg_arr


def guess_chromosome_regions(targets, telomere_size):
    """Determine (minimum) chromosome lengths from target coordinates."""
    endpoints = [subarr.end.iat[-1] for _c, subarr in targets.by_chromosome()]
    whole_chroms = GA.from_columns({
        'chromosome': targets.chromosome.drop_duplicates(),
        'start': telomere_size,
        'end': endpoints})
    return whole_chroms
