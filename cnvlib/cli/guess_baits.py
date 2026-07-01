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
single sample may have poor coverage at some targets by chance.
Still, this script does not guarantee correctly detecting all targets.

See also: https://github.com/brentp/goleft
"""

import argparse
import logging
import subprocess
import sys
from collections.abc import Iterable, Iterator
from typing import NamedTuple

import numpy as np
import pandas as pd

from skgenome import GenomicArray as GA
from skgenome import tabio

from .. import parallel
from ..coverage import do_coverage
from ..descriptives import modal_location


class Region(NamedTuple):
    """A captured genomic interval, 0-indexed half-open, with mean depth."""

    chromosome: str
    start: int
    end: int
    depth: float


def argument_parsing() -> argparse.Namespace:
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument(
        "sample_bams",
        nargs="+",
        help="""Sample BAM file(s) to test for target coverage.""",
    )
    AP.add_argument(
        "-o",
        "--output",
        metavar="FILENAME",
        help="""The inferred targets, in BED format.""",
    )
    AP.add_argument(
        "-c",
        "--coverage",
        metavar="FILENAME",
        help="""Filename to output average coverage depths in .cnn
                    format.""",
    )
    AP.add_argument(
        "-p",
        "--processes",
        metavar="CPU",
        nargs="?",
        type=int,
        const=0,
        default=1,
        help="""Number of subprocesses to segment in parallel.
                    If given without an argument, use the maximum number
                    of available CPUs. [Default: use 1 process]""",
    )
    AP.add_argument(
        "-f",
        "--fasta",
        metavar="FILENAME",
        help="Reference genome, FASTA format (e.g. UCSC hg19.fa)",
    )
    AP.add_argument(
        "-d",
        "--min-depth",
        metavar="DEPTH",
        type=int,
        default=5,
        help="""Minimum sequencing read depth to accept as captured. With
                    --access, half this depth is used for the initial scan.
                    [Default: %(default)s]""",
    )

    AP_x = AP.add_mutually_exclusive_group(required=True)
    AP_x.add_argument(
        "-t",
        "--targets",
        metavar="TARGET_BED",
        help="""Potentially targeted genomic regions, e.g. all known
                    exons in the reference genome, in BED format. Each of these
                    regions will be tested as a whole for enrichment. (Faster
                    method)""",
    )
    AP_x.add_argument(
        "-a",
        "--access",
        metavar="ACCESS_BED",
        help="""Sequencing-accessible genomic regions (e.g. from
                    'cnvkit.py access'), or known genic regions in the reference
                    genome, in BED format. All bases will be tested for
                    enrichment. (Slower method)""",
    )

    AP_access = AP.add_argument_group("With --access only")
    AP_access.add_argument(
        "-g",
        "--min-gap",
        metavar="GAP_SIZE",
        type=int,
        default=25,
        help="""Merge regions separated by gaps smaller than this.
                    [Default: %(default)s]""",
    )
    AP_access.add_argument(
        "-l",
        "--min-length",
        metavar="TARGET_SIZE",
        type=int,
        default=50,
        help="""Minimum region length to accept as captured.
                    [Default: %(default)s]""",
    )

    return AP.parse_args()


# ___________________________________________
# Guided method: guess from potential targets


def filter_targets(
    target_bed: str, sample_bams: list[str], procs: int, fasta: str | None
) -> GA:
    """Check if each potential target has significant coverage."""
    try:
        baits = tabio.read(target_bed, "bed4")
    except Exception as err:
        raise RuntimeError("Targets must be in BED format; try skg_convert.py") from err
    logging.info("Loaded %d candidate regions from %s", len(baits), target_bed)
    # Loop over BAMs to calculate weighted averages of bin coverage depths
    total_depths = np.zeros(len(baits), dtype=np.float64)
    for bam_fname in sample_bams:
        logging.info("Evaluating targets in %s", bam_fname)
        sample = do_coverage(target_bed, bam_fname, processes=procs, fasta=fasta)
        if len(sample) != len(baits):
            raise ValueError(
                f"Coverage returned {len(sample)} rows for {len(baits)} candidate "
                f"regions in {bam_fname}"
            )
        total_depths += sample["depth"].to_numpy()
    baits["depth"] = total_depths / len(sample_bams)
    logging.info("Average candidate-target depth:\n%s", baits["depth"].describe())
    return baits


# _________________________________________
# Unguided method: guess from raw depths


def scan_targets(
    access_bed: str,
    sample_bams: list[str],
    min_depth: float,
    min_gap: int,
    min_length: int,
    procs: int,
) -> GA:
    """Estimate baited regions from a genome-wide, per-base depth profile."""
    bait_chunks = []
    # ENH: context manager to call rm on bed chunks? with to_chunks as pool, ck?
    logging.info("Scanning for enriched regions in:\n  %s", "\n  ".join(sample_bams))
    with parallel.pick_pool(procs) as pool:
        args_iter = (
            (bed_chunk, sample_bams, min_depth, min_gap, min_length)
            for bed_chunk in parallel.to_chunks(access_bed)
        )
        for bed_chunk_fname, bait_chunk in pool.map(_scan_depth, args_iter):
            bait_chunks.append(bait_chunk)
            parallel.rm(bed_chunk_fname)
    baits = GA(pd.concat(bait_chunks))
    baits["depth"] /= len(sample_bams)
    return baits


def _scan_depth(args) -> tuple[str, pd.DataFrame]:
    """Wrapper for parallel map."""
    bed_fname, bam_fnames, min_depth, min_gap, min_length = args
    regions = list(
        drop_small(
            merge_gaps(scan_depth(bed_fname, bam_fnames, min_depth), min_gap),
            min_length,
        )
    )
    result = pd.DataFrame.from_records(regions, columns=Region._fields)
    if result.empty:
        # from_records infers object dtype for an empty frame, which would
        # down-cast every column to object in the scan_targets concat (most
        # genome-wide chunks are empty). Pin Region's dtypes so the numeric
        # columns survive and the --coverage log2 path stays valid.
        result = result.astype(
            {"chromosome": str, "start": int, "end": int, "depth": float}
        )
    return bed_fname, result


def scan_depth(
    bed_fname: str, bam_fnames: list[str], min_depth: float
) -> Iterator[Region]:
    """Locate sub-regions with enriched read depth in the given regions.

    Runs ``samtools depth`` over the BED regions and emits maximal runs of
    consecutive bases whose depth is at least ``min_depth``. A run is broken by
    a sub-threshold base, a chromosome change, or a gap in the reported
    positions: ``samtools depth`` omits bases with no confidently-mapped
    coverage -- uncaptured sequence, or the multi-mappers ``-Q 1`` drops (e.g.
    pseudogenes) -- which do not belong in the inferred panel. Each run is
    trimmed to the span where depth is at least half its peak.

    With multiple BAMs the per-base depths are summed across samples and
    ``min_depth`` is scaled by the number of samples.

    Yields
    ------
    Region
        A ``Region(chromosome, start, end, depth)`` namedtuple; coordinates are
        0-indexed half-open and ``depth`` is the mean over the trimmed span.
    """
    # NB: samtools emits additional BAMs' depths as trailing columns
    min_depth *= len(bam_fnames)

    with subprocess.Popen(
        [
            "samtools",
            "depth",
            "-Q",
            "1",  # Drop MAPQ-0 multi-mappers, keeping pseudogenes out of the panel
            "-b",
            bed_fname,
            *bam_fnames,
        ],
        stdout=subprocess.PIPE,
        text=True,
    ) as proc:
        # Detect runs of >= min_depth over consecutive positions; emit them
        chrom = ""
        start = 0
        prev_pos = 0
        depths: list[int] = []
        for line in proc.stdout:  # type: ignore[union-attr]
            fields = line.split("\t")
            pos = int(fields[1])
            depth = sum(map(int, fields[2:]))
            enriched = depth >= min_depth
            contiguous = fields[0] == chrom and pos == prev_pos + 1
            if depths and not (enriched and contiguous):
                # The current run ends at the previous base
                yield _make_region(chrom, start, depths)
                depths = []
            if enriched:
                if not depths:
                    # Open a new run
                    chrom = fields[0]
                    start = pos - 1
                depths.append(depth)
                prev_pos = pos
        # Flush the final run at end of stream
        if depths:
            yield _make_region(chrom, start, depths)

    if proc.returncode:
        raise RuntimeError(
            f"samtools depth failed (exit {proc.returncode}) scanning {bed_fname}"
        )


def _make_region(chrom: str, start: int, depths: list[int]) -> Region:
    """Trim a contiguous depth run to its >= half-peak span and average it."""
    darr = np.array(depths)
    half_depth = 0.5 * darr.max()
    ok_dp_idx = np.nonzero(darr >= half_depth)[0]
    start_idx = int(ok_dp_idx[0])
    end_idx = int(ok_dp_idx[-1]) + 1
    return Region(
        chrom,
        start + start_idx,
        start + end_idx,
        float(darr[start_idx:end_idx].mean()),
    )


def merge_gaps(regions: Iterable[Region], min_gap: int) -> Iterator[Region]:
    """Merge same-chromosome regions separated by gaps smaller than min_gap."""
    group: list[Region] = []
    for reg in regions:
        if group and (
            reg.chromosome != group[-1].chromosome
            or not 0 <= reg.start - group[-1].end < min_gap
        ):
            yield _merge_group(group)
            group = []
        group.append(reg)
    if group:
        yield _merge_group(group)


def _merge_group(group: list[Region]) -> Region:
    """Combine merged sub-regions; depth is their length-weighted mean."""
    if len(group) == 1:
        return group[0]
    first, last = group[0], group[-1]
    total_len = sum(reg.end - reg.start for reg in group)
    weighted_depth = sum(reg.depth * (reg.end - reg.start) for reg in group) / total_len
    return Region(first.chromosome, first.start, last.end, weighted_depth)


def drop_small(regions: Iterable[Region], min_length: int) -> Iterator[Region]:
    """Filter regions by minimum length."""
    return (reg for reg in regions if reg.end - reg.start >= min_length)


# ___________________________________________
# Shared


def normalize_depth_log2_filter(
    baits: GA, min_depth: float, enrich_ratio: float = 0.1
) -> GA:
    """Calculate normalized depth, add log2 column, filter by enrich_ratio."""
    # Normalize depths to a neutral value of 1.0
    dp_mode = modal_location(baits.data.loc[baits["depth"] > min_depth, "depth"].values)
    norm_depth = baits["depth"] / dp_mode
    # Drop low-coverage targets
    keep_idx = norm_depth >= enrich_ratio
    logging.info(
        "Keeping %d/%d bins with coverage depth >= %f, modal depth %f",
        keep_idx.sum(),
        len(keep_idx),
        dp_mode * enrich_ratio,
        dp_mode,
    )
    return baits[keep_idx]


def guess_baits(args) -> None:
    # Argparse yields 0 for "-p" with no argument; pick_pool reads <1 as
    # "use all available CPUs", so pass the process count through unchanged.
    if args.targets:
        baits = filter_targets(
            args.targets, args.sample_bams, args.processes, args.fasta
        )
    else:
        baits = scan_targets(
            args.access,
            args.sample_bams,
            0.5 * args.min_depth,  # More sensitive 1st pass
            args.min_gap,
            args.min_length,
            args.processes,
        )
    baits = normalize_depth_log2_filter(baits, args.min_depth)
    tabio.write(baits, args.output or sys.stdout, "bed")  # type: ignore[arg-type]
    if args.coverage:
        baits["log2"] = np.log2(baits["depth"] / baits["depth"].median())
        tabio.write(baits, args.coverage, "tab")


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    arguments = argument_parsing()
    guess_baits(arguments)


if __name__ == "__main__":
    main()
