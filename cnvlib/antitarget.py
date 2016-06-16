"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division, print_function
from builtins import map, next, range

import logging
import re
import sys

from . import tabio
from .params import INSERT_SIZE
from .gary import GenomicArray as GA


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
    target_chroms = dict(tabio.read_auto(target_bed).by_chromosome())
    if access_bed:
        # Chromosomes' accessible sequence regions are given -- use them
        access_chroms = dict(tabio.read_auto(access_bed).by_chromosome())
        if access_chroms and set(access_chroms).isdisjoint(target_chroms):
            raise ValueError("Chromosome names in the accessible regions file "
                             "%s %r do not match those in targets %s %r"
                             % (access_bed,
                                tuple(sorted(access_chroms.keys())[:3]),
                                target_bed,
                                tuple(sorted(target_chroms.keys())[:3])))
        # But filter out untargeted alternative contigs and mitochondria
        untgt_chroms = set(access_chroms) - set(target_chroms)
        # Autosomes typically have numeric names, allosomes are X and Y
        is_canonical = re.compile(r"(chr)?(\d+|[XYxy])$")
        if any(is_canonical.match(c) for c in target_chroms):
            chroms_to_skip = [c for c in untgt_chroms
                              if not is_canonical.match(c)]
        else:
            # Alternative contigs have long names -- skip them
            max_tgt_chr_name_len = max(list(map(len, target_chroms)))
            chroms_to_skip = [c for c in untgt_chroms
                              if len(c) > max_tgt_chr_name_len]
        for untgt_chr in chroms_to_skip:
            logging.info("Skipping untargeted chromosome %s", untgt_chr)
            del access_chroms[untgt_chr]
    else:
        # Chromosome accessible sequence regions not known -- use heuristics
        # (chromosome length is endpoint of last probe; skip initial
        # <magic number> of bases that are probably telomeric)
        TELOMERE_SIZE = 150000
        access_chroms = guess_chromosome_regions(target_chroms, TELOMERE_SIZE)

    backgrounds = find_background_regions(access_chroms, target_chroms,
                                          2 * INSERT_SIZE)
    bg_arr = GA.from_rows(backgrounds)
    bg_arr.sort()

    # Emit regions as antitarget bins according to avg_bin_size and min_bin_size
    out_rows = []
    for chrom, start, end in bg_arr.coords():
        span = end - start
        if span >= min_bin_size:
            nbins = int(round(span / avg_bin_size)) or 1
            if nbins == 1:
                out_rows.append((chrom, start, end))
            else:
                # Divide the background region into equal-sized bins
                bin_size = span / nbins
                bin_start = start
                bin_end = None
                for i in range(1, nbins):
                    bin_end = start + int(i * bin_size)
                    out_rows.append((chrom, bin_start, bin_end))
                    bin_start = bin_end
                out_rows.append((chrom, bin_start, end))
    out_arr = GA.from_rows(out_rows)
    out_arr["gene"] = "Background"
    return out_arr


def guess_chromosome_regions(target_chroms, telomere_size):
    """Determine (minimum) chromosome lengths from target coordinates."""
    endpoints = [target_region[len(target_region) - 1, 'end']
                 for _chrom, target_region in target_chroms.items()]
    whole_chroms = GA.from_columns({"chromosome": list(target_chroms.keys()),
                                    "start": telomere_size,
                                    "end": endpoints})
    return dict(whole_chroms.by_chromosome())


def find_background_regions(access_chroms, target_chroms, pad_size):
    """Take coordinates of accessible regions and targets; emit antitargets."""
    for chrom, access_arr in access_chroms.items():
        if chrom in target_chroms:
            target_regions = target_chroms[chrom].coords()

            # Split each access_region at the targets it contains
            _, tgt_start, tgt_end = next(target_regions)
            assert tgt_start < tgt_end
            for _, acc_start, acc_end in access_arr.coords():
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
                        _, tgt_start, tgt_end = next(target_regions)
                    except StopIteration:
                        # Ensure all the remaining access_regions are unbroken
                        tgt_start = tgt_end = float('Inf')

                # No remaining targets in this region - emit the rest of it
                if acc_end - pad_size - bg_start > 0:
                    yield (chrom, bg_start, acc_end - pad_size)
        else:
            for _, acc_start, acc_end in access_arr.coords():
                yield (chrom, acc_start + pad_size, acc_end - pad_size)

