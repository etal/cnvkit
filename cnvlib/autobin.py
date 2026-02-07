"""Estimate reasonable bin sizes from BAM read counts or depths."""

from __future__ import annotations

import logging
import os
import tempfile
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

import numpy as np
import pandas as pd
from skgenome import tabio, GenomicArray as GA

from . import coverage, samutil
from .antitarget import compare_chrom_names
from .descriptives import weighted_median


def midsize_file(fnames: Sequence[str]) -> str:
    """Select the median-size file from several given filenames.

    If an even number of files is given, selects the file just below the median.
    """
    assert fnames, "No files provided to calculate the median size."
    return sorted(fnames, key=lambda f: os.stat(f).st_size)[(len(fnames) - 1) // 2]


def do_autobin(
    bam_fname: str,
    method: str,
    targets: GA | None = None,
    access: GA | None = None,
    bp_per_bin: float = 100000.0,
    target_min_size: int = 20,
    target_max_size: int = 50000,
    antitarget_min_size: int = 500,
    antitarget_max_size: int = 1000000,
    fasta: str | None = None,
) -> tuple[tuple[float, int | None], tuple[float | None, int | None]]:
    """Quickly calculate reasonable bin sizes from BAM read counts.

    Parameters
    ----------
    bam_fname : str
        Path to BAM file.
    method : str
        Sequencing method: 'wgs' (whole-genome sequencing), 'amplicon' (targeted
        amplicon capture), or 'hybrid' (hybridization capture).
    targets : GenomicArray, optional
        Targeted genomic regions (required for 'hybrid' and 'amplicon').
    access : GenomicArray, optional
        Sequencing-accessible regions of the reference genome (for 'hybrid' and 'wgs').
    bp_per_bin : float, optional
        Desired number of sequencing read nucleotide bases mapped to each bin.
        Default is 100000.0.
    target_min_size : int, optional
        Minimum target bin size. Default is 20.
    target_max_size : int, optional
        Maximum target bin size. Default is 50000.
    antitarget_min_size : int, optional
        Minimum antitarget bin size. Default is 500.
    antitarget_max_size : int, optional
        Maximum antitarget bin size. Default is 1000000.
    fasta : str, optional
        Path to reference genome FASTA file.

    Returns
    -------
    tuple of tuple
        Nested tuple: ((target depth, target avg. bin size),
        (antitarget depth, antitarget avg. bin size)).

    Raises
    ------
    ValueError
        If targets are required for the method but not provided or empty.
    """
    if method in ("amplicon", "hybrid"):
        if targets is None:
            raise ValueError(
                f"Target regions are required for method {method!r} but were "
                "not provided."
            )
        if not len(targets):
            raise ValueError(
                f"Target regions are required for method {method!r} but were "
                "not provided."
            )

    # Closes over bp_per_bin
    def depth2binsize(depth: float | None, min_size: int, max_size: int) -> int | None:
        if not depth:
            return None
        bin_size = round(bp_per_bin / depth)
        if bin_size < min_size:
            logging.info(
                "Limiting est. bin size %d to given min. %d", bin_size, min_size
            )
            bin_size = min_size
        elif bin_size > max_size:
            logging.info(
                "Limiting est. bin size %d to given max. %d", bin_size, max_size
            )
            bin_size = max_size
        return bin_size

    samutil.ensure_bam_index(bam_fname)
    rc_table = samutil.idxstats(bam_fname, drop_unmapped=True, fasta=fasta)
    read_len = samutil.get_read_length(bam_fname, fasta=fasta)
    logging.info("Estimated read length %s", read_len)

    # Dispatch
    match method:
        case "amplicon":
            # From BAM index
            # rc_table = update_chrom_length(rc_table, targets)
            # tgt_depth = average_depth(rc_table, read_len)
            # By sampling
            assert targets is not None
            tgt_depth = sample_region_cov(bam_fname, targets, fasta=fasta)
            anti_depth: float | None = None
        case "hybrid":
            assert targets is not None
            tgt_depth, anti_depth = hybrid(
                rc_table,
                read_len,  # type: ignore[arg-type]
                bam_fname,
                targets,
                access,
                fasta,
            )
        case "wgs":
            if access is not None and len(access):
                rc_table = update_chrom_length(rc_table, access)
            tgt_depth = average_depth(rc_table, read_len)  # type: ignore[arg-type]
            anti_depth = None  # type: ignore[assignment]

    # Clip bin sizes to specified ranges
    tgt_bin_size = depth2binsize(tgt_depth, target_min_size, target_max_size)
    anti_bin_size = depth2binsize(anti_depth, antitarget_min_size, antitarget_max_size)
    return ((tgt_depth, tgt_bin_size), (anti_depth, anti_bin_size))


def hybrid(
    rc_table: pd.DataFrame,
    read_len: int | float,
    bam_fname: str,
    targets: GA,
    access: GA | None = None,
    fasta: str | None = None,
) -> tuple:
    """Hybrid capture sequencing."""
    # Identify off-target regions
    if access is None:
        access = idxstats2ga(rc_table, bam_fname)
        # Verify BAM chromosome names match those in target BED
        compare_chrom_names(access, targets)
    antitargets = access.subtract(targets)
    # Only examine chromosomes present in all 2-3 input datasets
    rc_table, targets, antitargets = shared_chroms(rc_table, targets, antitargets)
    # Deal with targets
    target_depth = sample_region_cov(bam_fname, targets, fasta=fasta)
    # Antitargets: subtract captured reads from total
    target_length = region_size_by_chrom(targets)["length"]
    target_reads = (target_length * target_depth / read_len).to_numpy()
    anti_table = update_chrom_length(rc_table, antitargets)
    anti_table = anti_table.assign(mapped=anti_table.mapped - target_reads)
    anti_depth = average_depth(anti_table, read_len)
    return target_depth, anti_depth


# ---


def average_depth(rc_table: pd.DataFrame, read_length: int | float) -> float:
    """Estimate the average read depth across the genome.

    Returns
    -------
    float
        Median of the per-chromosome mean read depths, weighted by chromosome
        size.
    """
    mean_depths = read_length * rc_table.mapped / rc_table.length
    return weighted_median(mean_depths, rc_table.length)  # type: ignore[no-any-return]


def idxstats2ga(table: pd.DataFrame, bam_fname: str) -> GA:
    return GA(
        table.assign(start=0, end=table.length).loc[:, ("chromosome", "start", "end")],
        meta_dict={"filename": bam_fname},
    )


def sample_region_cov(
    bam_fname: str, regions: GA, max_num: int = 100, fasta: str | None = None
) -> float:
    """Calculate read depth in a randomly sampled subset of regions."""
    midsize_regions = sample_midsize_regions(regions, max_num)
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w+t") as f:
        tabio.write(regions.as_dataframe(midsize_regions), f, "bed4")
        f.flush()
        table = coverage.bedcov(f.name, bam_fname, 0, fasta)
    # Mean read depth across all sampled regions
    return table.basecount.sum() / (table.end - table.start).sum()  # type: ignore[no-any-return]


def sample_midsize_regions(regions: GA, max_num: int) -> pd.DataFrame:
    """Randomly select a subset of up to `max_num` regions."""
    sizes = regions.end - regions.start
    lo_size, hi_size = np.percentile(sizes[sizes > 0], [25, 75])
    midsize_regions = regions.data[(sizes >= lo_size) & (sizes <= hi_size)]
    if len(midsize_regions) > max_num:
        midsize_regions = midsize_regions.sample(max_num, random_state=0xA5EED)
    return midsize_regions


def shared_chroms(*tables) -> list:
    """Intersection of DataFrame .chromosome values."""
    chroms = tables[0].chromosome.drop_duplicates()
    for tab in tables[1:]:
        if tab is not None:
            new_chroms = tab.chromosome.drop_duplicates()
            chroms = chroms[chroms.isin(new_chroms)]
    return [None if tab is None else tab[tab.chromosome.isin(chroms)] for tab in tables]


def update_chrom_length(rc_table: pd.DataFrame, regions: GA | None) -> pd.DataFrame:
    if regions is not None and len(regions):
        chrom_sizes = region_size_by_chrom(regions)
        rc_table = rc_table.merge(chrom_sizes, on="chromosome", how="inner")
        rc_table["length"] = rc_table["length_y"]  # ?
        rc_table = rc_table.drop(["length_x", "length_y"], axis=1)
    return rc_table


def region_size_by_chrom(regions: GA) -> pd.DataFrame:
    chromgroups = regions.data.groupby("chromosome", sort=False)
    # sizes = chromgroups.apply(total_region_size) # XXX
    sizes = [total_region_size(g) for _key, g in chromgroups]
    return pd.DataFrame(
        {"chromosome": regions.chromosome.drop_duplicates(), "length": sizes}
    )


def total_region_size(regions: GA) -> int:
    """Aggregate area of all genomic ranges in `regions`."""
    return (regions.end - regions.start).sum()  # type: ignore[no-any-return]
