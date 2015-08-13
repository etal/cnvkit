"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division
import sys
from itertools import groupby

from Bio._py3k import range
iteritems = (dict.iteritems if sys.version_info[0] < 3 else dict.items)

from . import core
from .params import INSERT_SIZE
from .rary import RegionArray as RA


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
    target_chroms = group_coords(RA.read(target_bed).coords())
    if access_bed:
        # Chromosome accessible sequence regions are given -- use them
        access_chroms = group_coords(RA.read(access_bed).coords())
    else:
        # Chromosome accessible sequence regions not known -- use heuristics
        # (chromosome length is endpoint of last probe; skip initial
        # <magic number> of bases that are probably telomeric)
        TELOMERE_SIZE = 150000
        access_chroms = guess_chromosome_regions(target_chroms, TELOMERE_SIZE)

    backgrounds = find_background_regions(access_chroms, target_chroms,
                                          2 * INSERT_SIZE)
    # Emit regions as antitarget bins according to avg_bin_size and min_bin_size
    # Do a set operation on backgrounds to avoid any duplicate regions
    for chrom, start, end in sorted(backgrounds, key=core.sorter_chrom_at(0)):
        span = end - start
        if span >= min_bin_size:
            nbins = int(round(span / avg_bin_size)) or 1
            if nbins == 1:
                yield (chrom, start, end)
            else:
                # Divide the background region into equal-sized bins
                bin_size = span / nbins
                bin_start = start
                bin_end = None
                for i in range(1, nbins):
                    bin_end = start + int(i * bin_size)
                    yield (chrom, bin_start, bin_end)
                    bin_start = bin_end
                yield (chrom, bin_start, end)


def group_coords(coordinates):
    """Group chromosomal coordinates into a dictionary."""
    out = {}
    for chrom, coords in groupby(coordinates, lambda cse: cse[0]):
        out[chrom] = [(start, end) for _chrom, start, end in coords]
    return out


def guess_chromosome_regions(target_chroms, telomere_size):
    """Determine (minimum) chromosome lengths from target coordinates."""
    out = {}
    for chrom, target_region in iteritems(target_chroms):
        endpoint = max(end for _start, end in target_region)
        # Put the tuple inside a list for compatibility w/ group_coords
        out[chrom] = [(telomere_size, endpoint)]
    return out


def find_background_regions(access_chroms, target_chroms, pad_size):
    """Take coordinates of accessible regions and targets; emit antitargets."""
    for chrom in access_chroms:
        if chrom in target_chroms:
            target_regions = iter(target_chroms[chrom])

            # Split each access_region at the targets it contains
            tgt_start, tgt_end = next(target_regions)
            assert tgt_start < tgt_end
            for acc_start, acc_end in access_chroms[chrom]:
                if acc_end - acc_start <= 2 * pad_size:
                    # Accessible region is way too small
                    continue

                bg_start = acc_start + pad_size
                while tgt_start < acc_end:
                    # There may be at least one more target in this region
                    if tgt_end + pad_size > bg_start:
                        # Yes, there is a target in this region
                        if tgt_start - pad_size > bg_start:
                            # Split the background region at this target
                            yield (chrom, bg_start, tgt_start - pad_size)
                        bg_start = tgt_end + pad_size

                    # Done splitting that target; is there another?
                    try:
                        tgt_start, tgt_end = next(target_regions)
                    except StopIteration:
                        # Ensure all the remaining access_regions are unbroken
                        tgt_start = tgt_end = float('Inf')

                # No remaining targets in this region - emit the rest of it
                if acc_end - pad_size - bg_start > 0:
                    yield (chrom, bg_start, acc_end - pad_size)
        else:
            for acc_start, acc_end in access_chroms[chrom]:
                yield (chrom, acc_start + pad_size, acc_end - pad_size)

