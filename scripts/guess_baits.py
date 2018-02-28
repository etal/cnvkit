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

import numpy as np
import pandas as pd

import cnvlib
from cnvlib import parallel
from cnvlib.descriptives import modal_location
from skgenome import tabio, GenomicArray as GA

logging.basicConfig(level=logging.INFO, format="%(message)s")


# ___________________________________________
# Guided method: guess from potential targets

def filter_targets(target_bed, sample_bams, procs):
    """Check if each potential target has significant coverage."""
    try:
        baits = tabio.read(target_bed, 'bed4')
    except:
        raise RuntimeError("Targets must be in BED format; try skg_convert.py")
    logging.info("Loaded %d candidate regions from %s", len(baits), target_bed)
    # Loop over BAMs to calculate weighted averages of bin coverage depths
    total_depths = np.zeros(len(baits), dtype=np.float_)
    for bam_fname in sample_bams:
        logging.info("Evaluating targets in %s", bam_fname)
        sample = cnvlib.do_coverage(target_bed, bam_fname, processes=procs)
        assert len(sample) == len(baits), \
                "%d != %d" % (len(sample), len(baits))
        total_depths += sample['depth'].values
    baits['depth'] = total_depths / len(sample_bams)
    logging.info("Average candidate-target depth:\n%s",
                 baits['depth'].describe())
    return baits


# _________________________________________
# Unguided method: guess from raw depths

def scan_targets(access_bed, sample_bams, min_depth, min_gap, min_length,
                 procs):
    """Estimate baited regions from a genome-wide, per-base depth profile."""
    bait_chunks = []
    # ENH: context manager to call rm on bed chunks? with to_chunks as pool, ck?
    logging.info("Scanning for enriched regions in:\n  %s",
                 '\n  '.join(sample_bams))
    #  with futures.ProcessPoolExecutor(procs) as pool:
    with parallel.pick_pool(procs) as pool:
        args_iter = ((bed_chunk, sample_bams,
                    min_depth, min_gap, min_length)
                    for bed_chunk in parallel.to_chunks(access_bed))
        for bed_chunk_fname, bait_chunk in pool.map(_scan_depth, args_iter):
            bait_chunks.append(bait_chunk)
            parallel.rm(bed_chunk_fname)
    baits = GA(pd.concat(bait_chunks))
    baits['depth'] /= len(sample_bams)
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

    nsamples = len(bam_fnames)
    if nsamples == 1:
        def get_depth(depths):
            return int(depths[0])
    else:
        min_depth *= nsamples
        # NB: samtools emits additional BAMs' depths as trailing columns
        def get_depth(depths):
            return sum(map(int, depths))

    proc = subprocess.Popen(['samtools', 'depth',
                             '-Q',  '1',  # Skip pseudogenes
                             '-b', bed_fname,
                            ] + bam_fnames,
                            stdout=subprocess.PIPE,
                            shell=False)

    # Detect runs of >= min_depth; emit their coordinates
    chrom = start = depths = None
    for line in proc.stdout:
        fields = line.split('\t')
        depth = get_depth(fields[2:])
        is_enriched = (depth >= min_depth)
        if start is None:
            if is_enriched:
                # Entering a new captured region
                chrom = fields[0]
                start = int(fields[1]) - 1
                depths = [depth]
            else:
                # Still not in a captured region
                continue
        elif is_enriched and fields[0] == chrom:
            # Still in a captured region -- extend it
            depths.append(depth)
        else:
            # Exiting a captured region
            # Update target region boundaries
            darr = np.array(depths)
            half_depth = 0.5 * darr.max()
            ok_dp_idx = np.nonzero(darr >= half_depth)[0]
            start_idx = ok_dp_idx[0]
            end_idx = ok_dp_idx[-1] + 1
            yield Region(chrom,
                            start + start_idx,
                            start + end_idx,
                            darr[start_idx:end_idx].mean())
            chrom = start = depths = None


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


# ___________________________________________
# Shared

def normalize_depth_log2_filter(baits, min_depth, enrich_ratio=0.1):
    """Calculate normalized depth, add log2 column, filter by enrich_ratio."""
    # Normalize depths to a neutral value of 1.0
    dp_mode = modal_location(baits.data.loc[baits['depth'] > min_depth,
                                            'depth'].values)
    norm_depth = baits['depth'] / dp_mode
    # Drop low-coverage targets
    keep_idx = (norm_depth >= enrich_ratio)
    logging.info("Keeping %d/%d bins with coverage depth >= %f, modal depth %f",
                 keep_idx.sum(), len(keep_idx), dp_mode * enrich_ratio, dp_mode)
    return baits[keep_idx]



if __name__ == '__main__':
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('sample_bams', nargs='+',
                    help="""Sample BAM file(s) to test for target coverage.""")
    AP.add_argument('-o', '--output', metavar='FILENAME',
                    help="""The inferred targets, in BED format.""")
    AP.add_argument('-c', '--coverage', metavar='FILENAME',
                    help="""Filename to output average coverage depths in .cnn
                    format.""")
    AP.add_argument('-p', '--processes', metavar='CPU',
                    nargs='?', type=int, const=0, default=1,
                    help="""Number of subprocesses to segment in parallel.
                    If given without an argument, use the maximum number
                    of available CPUs. [Default: use 1 process]""")

    AP_x = AP.add_mutually_exclusive_group(required=True)
    AP_x.add_argument('-t', '--targets', metavar='TARGET_BED',
                    help="""Potentially targeted genomic regions, e.g. all known
                    exons in the reference genome, in BED format. Each of these
                    regions will be tested as a whole for enrichment. (Faster
                    method)""")
    AP_x.add_argument('-a', '--access', metavar='ACCESS_BED',
                    # default="../data/access-5k-mappable.grch37.bed",
                    help="""Sequencing-accessible genomic regions (e.g. from
                    'cnvkit.py access'), or known genic regions in the reference
                    genome, in BED format. All bases will be tested for
                    enrichment. (Slower method)""")

    AP_target = AP.add_argument_group("With --targets only")
    AP_target.add_argument('-d', '--min-depth', metavar='DEPTH',
                    type=int, default=5,
                    help="""Minimum sequencing read depth to accept as captured.
                    [Default: %(default)s]""")

    AP_access = AP.add_argument_group("With --access only")
    AP_access.add_argument('-g', '--min-gap', metavar='GAP_SIZE',
                    type=int, default=25,
                    help="""Merge regions separated by gaps smaller than this.
                    [Default: %(default)s]""")
    AP_access.add_argument('-l', '--min-length', metavar='TARGET_SIZE',
                    type=int, default=50,
                    help="""Minimum region length to accept as captured.
                    [Default: %(default)s]""")

    args = AP.parse_args()

    # ENH: can we reserve multiple cores for htslib?
    if args.processes < 1:
        args.processes = None

    if args.targets:
        baits = filter_targets(args.targets, args.sample_bams, args.processes)
    else:
        baits = scan_targets(args.access, args.sample_bams,
                             0.5 * args.min_depth,  # More sensitive 1st pass
                             args.min_gap, args.min_length, args.processes)
    baits = normalize_depth_log2_filter(baits, args.min_depth)
    tabio.write(baits, args.output or sys.stdout, 'bed')
    if args.coverage:
        baits['log2'] = np.log2(baits['depth'] / baits['depth'].median())
        tabio.write(baits, args.coverage, 'tab')
