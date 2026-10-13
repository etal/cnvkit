"""Estimate reasonable bin sizes from read counts or depths.

The sample may be a BAM or CRAM file, or a tabix-indexed bedGraph (.bed.gz) of
per-base read depth.
"""

from __future__ import annotations

import logging
import os
import tempfile
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

import numpy as np
import pandas as pd

from skgenome import GenomicArray as GA
from skgenome import tabio

from . import coverage, samutil
from .antitarget import compare_chrom_names
from .descriptives import weighted_median

# A bedGraph is read line by line, so each region sampled from it is trimmed to a
# centered window of at most this many bases. Accessible and antitarget regions
# span up to megabases each, and 100 of them untrimmed would mean reading tens of
# megabases of per-base depth; targets are almost always shorter than this.
BEDGRAPH_SAMPLE_WIDTH = 10000


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
    """Quickly calculate reasonable bin sizes from read counts or depths.

    Parameters
    ----------
    bam_fname : str
        Path to a BAM or CRAM file, or to a tabix-indexed bedGraph file
        (.bed.gz) of per-base read depth.
    method : str
        Sequencing method: 'wgs' (whole-genome sequencing), 'amplicon' (targeted
        amplicon capture), or 'hybrid' (hybridization capture).
    targets : GenomicArray, optional
        Targeted genomic regions (required for 'hybrid' and 'amplicon').
    access : GenomicArray, optional
        Sequencing-accessible regions of the reference genome (for 'hybrid' and
        'wgs'). Required for 'wgs' on a bedGraph, and for a bedGraph's 'hybrid'
        antitarget depth.
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
        If targets are required for the method but not provided or empty; if
        method 'wgs' is given a bedGraph without accessible regions; or if no
        read depth is found in the target regions.

    Notes
    -----
    A BAM or CRAM file is summarized by its index: per-chromosome read counts
    and the estimated read length give the 'wgs' depth and the 'hybrid'
    antitarget depth. A bedGraph records neither read counts nor contig
    lengths, so every depth is instead measured in a seeded random sample of
    regions -- the targets, the accessible regions ('wgs'), or the accessible
    regions on targeted chromosomes minus the targets ('hybrid' antitargets) --
    as the target depth is for any input. The 'hybrid' target depth of a
    bedGraph covers all targets, whereas for a BAM it covers only those on
    accessible chromosomes.
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
        if depth is None or not np.isfinite(depth) or depth <= 0:
            # Guards the antitarget path (anti_depth may be None, or NaN when no
            # reads fall off-target): a bare `if not depth` misses NaN, since
            # `not np.nan` is False, letting it reach round() -> ValueError.
            # (do_autobin already raises on a bad *target* depth before here.)
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

    if access is not None and not len(access):
        # Empty accessible regions mean none given, as on the command line
        access = None
    if coverage.is_bedgraph_format(bam_fname):
        tgt_depth, anti_depth = bedgraph_depths(bam_fname, method, targets, access)
    else:
        tgt_depth, anti_depth = bam_depths(bam_fname, method, targets, access, fasta)

    # A non-finite or non-positive target depth means no reads were found in
    # the sampled/accessible regions -- give an actionable error rather than a
    # cryptic "cannot convert float NaN to integer" downstream (#421).
    if tgt_depth is None or not np.isfinite(tgt_depth) or tgt_depth <= 0:
        raise ValueError(
            f"Could not estimate read depth in the target regions of "
            f"{bam_fname!r} (method {method!r}): got {tgt_depth!r}. This "
            "usually means no reads mapped to those regions -- most often "
            "because the chromosome names in the sample file don't match those "
            "in the reference, --access, or --targets regions (e.g. 'chr1' "
            "vs. '1'), or the regions are empty. Check that the sample file, "
            "reference FASTA, and any --access/--targets files use consistent "
            "chromosome names."
        )

    # Clip bin sizes to specified ranges
    tgt_bin_size = depth2binsize(tgt_depth, target_min_size, target_max_size)
    anti_bin_size = depth2binsize(anti_depth, antitarget_min_size, antitarget_max_size)
    return ((tgt_depth, tgt_bin_size), (anti_depth, anti_bin_size))


def bam_depths(
    bam_fname: str,
    method: str,
    targets: GA | None,
    access: GA | None,
    fasta: str | None,
) -> tuple[float, float | None]:
    """Target and antitarget read depths of a BAM or CRAM file."""
    samutil.ensure_bam_index(bam_fname)
    rc_table = samutil.idxstats(bam_fname, drop_unmapped=True, fasta=fasta)
    read_len = samutil.get_read_length(bam_fname, fasta=fasta)
    logging.info("Estimated read length %s", read_len)

    # Dispatch
    match method:
        case "amplicon":
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
            rc_table = update_chrom_length(rc_table, access)
            tgt_depth = average_depth(rc_table, read_len)  # type: ignore[arg-type]
            anti_depth = None
        case _:
            raise ValueError(f"Unknown sequencing method {method!r}")
    return tgt_depth, anti_depth


def bedgraph_depths(
    bg_fname: str, method: str, targets: GA | None, access: GA | None
) -> tuple[float, float | None]:
    """Target and antitarget read depths of a bedGraph, by region sampling."""
    match method:
        case "amplicon":
            assert targets is not None
            return sample_region_cov(bg_fname, targets), None
        case "hybrid":
            assert targets is not None
            tgt_depth = sample_region_cov(bg_fname, targets)
            if access is None:
                logging.warning(
                    "Antitarget depth can't be estimated from bedGraph %s without "
                    "sequencing-accessible regions (--access); a bedGraph does "
                    "not record chromosome lengths",
                    bg_fname,
                )
                return tgt_depth, None
            # Off-target regions on the targeted chromosomes, as hybrid() takes
            # for BAM. Unlike hybrid(), which subtracts target reads from each
            # chromosome's total, the target depth needs no such restriction.
            antitargets = shared_chroms(access.subtract(targets), targets)[0]
            if not len(antitargets):
                logging.warning(
                    "No accessible regions lie off-target on the targeted "
                    "chromosomes; check that the --access and --targets "
                    "chromosome names match. Antitarget depth can't be estimated."
                )
                return tgt_depth, None
            return tgt_depth, sample_region_cov(bg_fname, antitargets)
        case "wgs":
            if access is None:
                raise ValueError(
                    f"Method 'wgs' needs sequencing-accessible regions (--access) "
                    f"to estimate the read depth of bedGraph {bg_fname!r}, since a "
                    "bedGraph does not record chromosome lengths."
                )
            return sample_region_cov(bg_fname, access), None
        case _:
            raise ValueError(f"Unknown sequencing method {method!r}")


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
        size. NaN if no chromosomes remain (e.g. the BAM and region chromosome
        names don't overlap), which the caller turns into an actionable error.
    """
    if rc_table.empty:
        return float("nan")
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
    """Calculate read depth in a randomly sampled subset of regions.

    `bam_fname` may be a BAM or CRAM file, or a bedGraph (.bed.gz), whose
    sampled regions are first trimmed to `BEDGRAPH_SAMPLE_WIDTH` bases.
    """
    if not (regions.end > regions.start).any():
        # No targets on the accessible chromosomes, or only zero-width regions:
        # depth is undefined, and do_autobin raises an actionable error.
        return float("nan")
    midsize_regions = sample_midsize_regions(regions, max_num)
    is_bedgraph = coverage.is_bedgraph_format(bam_fname)
    if is_bedgraph:
        midsize_regions = center_windows(midsize_regions, BEDGRAPH_SAMPLE_WIDTH)
    else:
        bam_chroms = samutil.get_bam_chroms(bam_fname, fasta)
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w+t") as f:
        tabio.write(regions.as_dataframe(midsize_regions), f, "bed4")
        f.flush()
        if is_bedgraph:
            table = coverage.bedgraph_to_basecount(f.name, bam_fname)
        else:
            table = coverage.bedcov(
                f.name, bam_fname, 0, fasta, allow_empty=True, bam_chroms=bam_chroms
            )
    if not len(table):
        # No sampled region overlaps a BAM contig -- depth is undefined; let
        # the caller (do_autobin) raise an actionable chrom-name-mismatch error.
        return float("nan")
    # Mean read depth across all sampled regions
    return table.basecount.sum() / (table.end - table.start).sum()  # type: ignore[no-any-return]


def sample_midsize_regions(regions: GA, max_num: int) -> pd.DataFrame:
    """Randomly select a subset of up to `max_num` regions.

    Regions between the 25th and 75th percentile of size are preferred. When no
    region falls in that range, as with two regions of very different size, all
    regions of nonzero size are eligible instead.
    """
    sizes = regions.end - regions.start
    lo_size, hi_size = np.percentile(sizes[sizes > 0], [25, 75])
    in_range = (sizes >= lo_size) & (sizes <= hi_size)
    midsize_regions = regions.data[in_range if in_range.any() else sizes > 0]
    if len(midsize_regions) > max_num:
        midsize_regions = midsize_regions.sample(max_num, random_state=0xA5EED)
    return midsize_regions


def center_windows(regions: pd.DataFrame, width: int) -> pd.DataFrame:
    """Trim each region longer than `width` to its centered `width`-base window."""
    lengths = regions.end - regions.start
    starts = regions.start + (lengths - width).clip(lower=0) // 2
    return regions.assign(start=starts, end=starts + lengths.clip(upper=width))


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
