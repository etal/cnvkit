#!/usr/bin/env python
"""Guess the coordinates of captured regions from sample read depths.

Two approaches available:

- (Faster) Scan a given list of exons and/or other potentially targeted regions.
  The script checks each region and drops those with very low coverage
  indicating they were not captured.
- (Slower) Scan the entire genome, or the given sequencing-accessible regions,
  for regions with elevated coverage. Choose reasonable boundaries for each
  apparently captured region.

Use multiple BAMs for greater robustness in detecting targeted regions, as a
single sample may have poor coverage are some targets by chance.
Still, this script does not guarantee correctly detecting all targets.

See also: https://github.com/brentp/goleft
"""
from __future__ import absolute_import, division, print_function

import argparse
import collections
import logging
import subprocess
import sys
from concurrent import futures

import numpy as np
import pandas as pd

import cnvlib
from cnvlib import tabio
from cnvlib.parallel import to_chunks, rm
from cnvlib.descriptives import modal_location
from cnvlib.genome import GenomicArray as GA

logging.basicConfig(level=logging.INFO, format="%(message)s")


# ___________________________________________
# Guided method: guess from potential targets

def filter_targets(target_bed, sample_bams, min_depth, procs):
    """Check if each potential target has significant coverage."""
    baits = tabio.read_auto(target_bed)
    # Loop over BAMs to calculate weighted averages of bin coverage depths
    total_depths = np.zeros(len(baits), dtype=np.float_)
    for bam_fname in sample_bams:
        logging.info("Scanning targets in %s", bam_fname)
        sample = cnvlib.do_coverage(target_bed, bam_fname, processes=procs)
        total_depths += sample['depth'].values
    min_total_depth = min_depth * len(sample_bams)
    # Normalize depths to a neutral value of 1.0
    mode_dp = modal_location(total_depths[total_depths > min_total_depth])
    baits['depth'] = total_depths / len(sample_bams)
    norm_depth = total_depths / mode_dp
    baits['log2'] = np.log2(norm_depth)
    # Drop low-coverage targets
    ok_idx = (norm_depth >= 0.1)  # Require at least 10x enrichment
    logging.info("Keeping %d/%d bins with average coverage depth >= %f, "
                 "(initially %f), mode %f",
                 ok_idx.sum(), len(ok_idx), mode_dp * 0.1, min_depth, mode_dp)
    return baits[ok_idx]


# _________________________________________
# Unguided method: guess from raw depths

def scan_targets(access_bed, sample_bams, min_depth, min_gap, min_length,
                 procs):
    """Estimate baited regions from a genome-wide, per-base depth profile."""
    bait_chunks = []
    # ENH: context manager to call rm on bed chunks?
    logging.info("Scanning for enriched regions in:\n  %s",
                 '\n  '.join(sample_bams))
    with futures.ProcessPoolExecutor(procs) as pool:
        args_iter = ((bed_chunk, sample_bams,
                    min_depth, min_gap, min_length)
                    for bed_chunk in to_chunks(access_bed))
        for bed_chunk_fname, bait_chunk in pool.map(_scan_depth, args_iter):
            bait_chunks.append(bait_chunk)
            rm(bed_chunk_fname)
    baits = GA(pd.concat(bait_chunks))
    return baits


def _scan_depth(args):
    """Wrapper for parallel map"""
    bed_fname, bam_fnames, min_depth, min_gap, min_length = args
    regions = list(drop_small(merge_gaps(scan_depth(bed_fname, bam_fnames,
                                                    min_depth),
                                         min_gap),
                              min_length))
    result = pd.DataFrame.from_records(list(regions),
                                       columns=regions[0]._fields)
    return bed_fname, result


def scan_depth(bed_fname, bam_fnames, min_depth):
    """Locate sub-regions with enriched read depth in the given regions.

    Yields
    ------
    tuple
        Region coordinates (0-indexed, half-open): chromosome name, start, end
    """
    Region = collections.namedtuple('Region', 'chromosome start end depth')

    if len(bam_fnames) == 1:
        def get_depth(depths):
            return int(depths[0])
    else:
        # ENH: take multiple BAMs; samtools emits add'l cols
        def get_depth(depths):
            return np.median(map(int, depths))

    proc = subprocess.Popen(['samtools', 'depth',
                             '-Q',  '1',  # Skip pseudogenes
                             '-b', bed_fname,
                            ] + bam_fnames,
                            stdout=subprocess.PIPE,
                            shell=False)

    # Detect runs of >= min_depth; emit their coordinates
    # ENH: Use the depth profile within each candidate region to choose a
    # cutoff for the bait boundary, e.g. 1/2 the region's peak.
    curr_chrom = start = end = None
    tot_depth = 0
    for line in proc.stdout:
        fields = line.split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        depth = get_depth(fields[2:])
        is_enriched = (depth >= min_depth)
        if start is None:
            if is_enriched:
                # Entering a new captured region
                curr_chrom = chrom
                start = pos - 1
                end = pos
                tot_depth = depth
            else:
                # Still not in a captured region
                continue
        else:
            if is_enriched and chrom == curr_chrom:
                # Still in a captured region -- extend it
                end = pos
                tot_depth += depth
            else:
                # Exiting a captured region
                yield Region(curr_chrom, start, end,
                             tot_depth / (end - start))
                curr_chrom = start = end = None
                tot_depth = 0


def merge_gaps(regions, min_gap):
    """Merge regions across small gaps."""
    prev = next(regions)
    for reg in regions:
        if reg.start - prev.end < min_gap:
            # Merge
            prev = prev._replace(end=reg.end)
        else:
            # Done merging; move on
            yield prev
            prev = reg
    # Residual
    yield prev


def drop_small(regions, min_length):
    """Merge small gaps and filter by minimum length."""
    return (reg for reg in regions
            if reg.end - reg.start >= min_length)



if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('sample_bams', nargs='+',
                    # default="R15-08280_kapa-NGv3-PE100-NGv3-ready.bam",
                    help="""Sample BAM file to test for target coverage""")
    AP.add_argument('-o', '--output')
    AP.add_argument('-c', '--coverage',
                    help="""Filename to output average coverage depths in .cnn
                    format.""")
    AP.add_argument('-p', '--processes', nargs='?', type=int, const=0, default=1,
                    help="""Number of subprocesses to segment in parallel.
                    If given without an argument, use the maximum number
                    of available CPUs. [Default: use 1 process]""")

    AP_x = AP.add_mutually_exclusive_group(required=True)
    AP_x.add_argument('-t', '--targets', metavar='TARGET_BED',
                    help="""Potentially targeted genomic regions, e.g. all
                    possible exons for the reference genome.""")
    AP_x.add_argument('-a', '--access', metavar='ACCESS_BED',
                    # default="~/code/cnvkit/data/access-5k-mappable.grch37.bed",
                    help="""Sequencing-accessible genomic regions, or exons to
                    use as possible targets (e.g. output of refFlat2bed.py)""")

    AP.add_argument('-d', '--min-depth', type=int, default=5,
                    help="Minimum sequencing read depth to accept as captured.")
    # Just for --access
    AP.add_argument('-g', '--min-gap', type=int, default=25,
                    help="Merge regions with gaps below this threshold.")
    AP.add_argument('-l', '--min-length', type=int, default=100,
                    help="Minimum region length to accept as captured.")
    # ENH: for -t, another option to include each intron in testing

    args = AP.parse_args()

    # ENH: can we reserve multiple cores for htslib?
    if args.processes < 1:
        args.processes = None

    if args.targets:
        baits = filter_targets(args.targets, args.sample_bams, args.min_depth,
                               args.processes)
    else:
        baits = scan_targets(args.access, args.sample_bams, args.min_depth,
                             args.min_gap, args.min_length, args.processes)
    tabio.write(baits, args.output or sys.stdout, 'bed')
    if args.coverage:
        baits['log2'] = np.log2(baits['depth'] / baits['depth'].median())
        tabio.write(baits, args.coverage, 'tab')
