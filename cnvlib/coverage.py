"""Supporting functions for the 'antitarget' command."""

from __future__ import annotations
import gzip
import logging
import math
import os.path
import tempfile
import time
from concurrent import futures
from io import StringIO
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from skgenome import tabio
from skgenome._pysam import PYSAM_INSTALL_MSG

from . import core, samutil
from .cnary import CopyNumArray as CNA
from .parallel import rm, to_chunks
from .params import NULL_LOG2_COVERAGE

if TYPE_CHECKING:
    import pysam
    from collections.abc import Iterator
    from skgenome import GenomicArray


# SAM flags for reads excluded from coverage: unmapped (0x4), secondary (0x100),
# QC-fail (0x200), and duplicate (0x400), i.e. 0x704. This mirrors the
# `filter_read` predicate used by the by-count path (`region_depth_count`) and
# happens to equal `samtools bedcov`'s own default exclude set -- but we pass it
# to bedcov explicitly so the depth path keeps ignoring marked duplicates even
# if a samtools version changes that default (#689). Picard MarkDuplicates sets
# the duplicate flag, so marking duplicates is equivalent to removing them.
COVERAGE_EXCLUDE_FLAGS = 0x704


def is_bedgraph_format(fname: str) -> bool:
    """Check if input file is bedGraph format (.bed.gz)."""
    return fname.endswith(".bed.gz")


def validate_bedgraph_index(fname: str) -> None:
    """Ensure bedGraph file has a tabix index.

    Raises
    ------
    FileNotFoundError
        If the index file (.tbi or .csi) is not found.
    """
    if not (os.path.exists(f"{fname}.tbi") or os.path.exists(f"{fname}.csi")):
        raise FileNotFoundError(
            f"bedGraph file {fname} requires a tabix index (.tbi or .csi). "
            f"Create one with: tabix -p bed {fname}"
        )


def bedgraph_to_basecount(
    bed_fname: str,
    bedgraph_fname: str,
) -> pd.DataFrame:
    """Calculate coverage in BED regions from a bedGraph file.

    Reads a tabix-indexed bedGraph file and calculates the total depth
    (sum of depth * bases) for each region in the BED file.

    Parameters
    ----------
    bed_fname : str
        Path to BED file defining regions to measure coverage.
    bedgraph_fname : str
        Path to tabix-indexed bedGraph file (.bed.gz).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns matching bedcov output: chromosome, start,
        end, gene (if present in BED), and basecount.
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            f"pysam is required for reading bedGraph files. {PYSAM_INSTALL_MSG}"
        ) from None

    validate_bedgraph_index(bedgraph_fname)

    # Read the BED file to get regions
    bed_data = tabio.read(bed_fname, "bed")

    # Open tabix file for random access with guaranteed cleanup
    try:
        tbx = pysam.TabixFile(bedgraph_fname)
    except OSError as exc:
        raise ValueError(
            f"Failed to open bedGraph file {bedgraph_fname!r}. Error: {exc}"
        ) from exc

    try:
        rows = []
        chromosomes_in_bedgraph = set(tbx.contigs)

        for region in bed_data:
            chrom = region.chromosome
            start = region.start
            end = region.end

            # Handle chromosome naming mismatches or missing chromosomes
            if chrom not in chromosomes_in_bedgraph:
                # Try adding/removing 'chr' prefix
                alt_chrom = (
                    f"chr{chrom}"
                    if not chrom.startswith("chr")
                    else chrom.removeprefix("chr")
                )
                if alt_chrom in chromosomes_in_bedgraph:
                    chrom = alt_chrom
                else:
                    # Chromosome not in bedGraph - treat as 0 coverage
                    logging.debug(
                        "Chromosome %s not found in bedGraph %s, using 0 coverage",
                        region.chromosome,
                        bedgraph_fname,
                    )
                    basecount = 0.0
                    row_data = {
                        "chromosome": region.chromosome,
                        "start": start,
                        "end": end,
                        "basecount": basecount,
                    }
                    if hasattr(region, "gene") and region.gene:
                        row_data["gene"] = region.gene
                    rows.append(row_data)
                    continue

            # Query bedGraph for this region
            basecount = 0.0
            try:
                for line in tbx.fetch(chrom, start, end):
                    fields = line.split("\t")
                    r_start = int(fields[1])
                    r_end = int(fields[2])
                    depth = float(fields[3])

                    # Clip to bin boundaries
                    overlap_start = max(r_start, start)
                    overlap_end = min(r_end, end)
                    if overlap_end > overlap_start:
                        basecount += depth * (overlap_end - overlap_start)
            except Exception as exc:
                logging.warning(
                    "Error querying bedGraph for %s:%d-%d: %s",
                    chrom,
                    start,
                    end,
                    exc,
                )

            # Build row matching bedcov output format
            row_data = {
                "chromosome": region.chromosome,  # Use original name from BED
                "start": start,
                "end": end,
                "basecount": basecount,
            }
            if hasattr(region, "gene") and region.gene:
                row_data["gene"] = region.gene
            rows.append(row_data)

        # Create DataFrame with columns matching bedcov output
        if rows and "gene" in rows[0]:
            columns = ["chromosome", "start", "end", "gene", "basecount"]
        else:
            columns = ["chromosome", "start", "end", "basecount"]

        return pd.DataFrame(rows, columns=columns)
    finally:
        tbx.close()


def do_coverage(
    bed_fname: str,
    bam_or_bg_fname: str,
    by_count: bool = False,
    min_mapq: int = 0,
    processes: int = 1,
    fasta: str | None = None,
) -> CNA:
    """Calculate coverage in the given regions from BAM read depths.

    Parameters
    ----------
    bed_fname : str
        Path to BED file defining regions to measure coverage.
    bam_or_bg_fname : str
        Path to BAM file containing aligned reads or bedGraph file (.bed.gz).
    by_count : bool, optional
        Calculate coverage by read count instead of read depth. Default is False.
        Ignored for bedGraph input.
    min_mapq : int, optional
        Minimum mapping quality score to include a read. Default is 0.
        Ignored for bedGraph input.
    processes : int, optional
        Number of parallel processes to use. Default is 1.
        Ignored for bedGraph input.
    fasta : str, optional
        Path to reference genome FASTA file.

    Returns
    -------
    CopyNumArray
        Coverage values for each region.

    Raises
    ------
    RuntimeError
        If the BAM file is not sorted by coordinates.
    FileNotFoundError
        If bedGraph file lacks required tabix index.
    """
    # Detect input format
    if is_bedgraph_format(bam_or_bg_fname):
        # bedGraph format - use simplified processing
        if by_count:
            logging.warning(
                "Option --count/-c is not applicable to bedGraph input and will be ignored"
            )
        if min_mapq > 0:
            logging.warning(
                "Option --min-mapq/-q is not applicable to bedGraph input and will be ignored"
            )
        if processes > 1:
            logging.warning(
                "Option --processes/-p is not applicable to bedGraph input and will be ignored"
            )
        cnarr = interval_coverages_bedgraph(bed_fname, bam_or_bg_fname)
    else:
        # BAM format - use existing logic
        if not samutil.ensure_bam_sorted(bam_or_bg_fname, fasta=fasta):
            raise RuntimeError(
                f"BAM file {bam_or_bg_fname} must be sorted by coordinates"
            )
        samutil.ensure_bam_index(bam_or_bg_fname)
        cnarr = interval_coverages(
            bed_fname, bam_or_bg_fname, by_count, min_mapq, processes, fasta
        )

    return cnarr


def interval_coverages(
    bed_fname: str,
    bam_fname: str,
    by_count: bool,
    min_mapq: int,
    processes: int,
    fasta: str | None = None,
) -> CNA:
    """Calculate log2 coverages in the BAM file at each interval."""
    meta = {"sample_id": core.fbase(bam_fname)}
    start_time = time.time()

    # Skip processing if the BED file is empty
    with open(bed_fname) as bed_handle:
        for line in bed_handle:
            if line.strip():
                break
        else:
            logging.info(
                "Skip processing %s with empty regions file %s",
                os.path.basename(bam_fname),
                bed_fname,
            )
            return CNA.from_rows([], meta_dict=meta)  # type: ignore[return-value]

    # Calculate average read depth in each bin
    if by_count:
        results = interval_coverages_count(
            bed_fname,
            bam_fname,
            min_mapq,
            processes,
            fasta,  # type: ignore[arg-type]
        )
        read_counts, cna_rows = zip(*results, strict=True)
        read_counts = pd.Series(read_counts)
        cnarr = CNA.from_rows(
            list(cna_rows), columns=(*CNA._required_columns, "depth"), meta_dict=meta
        )
    else:
        table = interval_coverages_pileup(
            bed_fname, bam_fname, min_mapq, processes, fasta
        )
        read_len = samutil.get_read_length(bam_fname, fasta=fasta)
        read_counts = table["basecount"] / read_len
        table = table.drop("basecount", axis=1)
        cnarr = CNA(table, meta)

    # Log some stats
    tot_time = time.time() - start_time
    tot_reads = read_counts.sum()
    logging.info(
        "Time: %.3f seconds (%d reads/sec, %s bins/sec)",
        tot_time,
        int(round(tot_reads / tot_time, 0)),
        int(round(len(read_counts) / tot_time, 0)),
    )
    logging.info(
        "Summary: #bins=%d, #reads=%d, mean=%.4f, min=%s, max=%s",
        len(read_counts),
        tot_reads,
        (tot_reads / len(read_counts)),
        read_counts.min(),
        read_counts.max(),
    )
    tot_mapped_reads = samutil.bam_total_reads(bam_fname, fasta=fasta)
    if tot_mapped_reads:
        logging.info(
            "Percent reads in regions: %.3f (of %d mapped)",
            100.0 * tot_reads / tot_mapped_reads,
            tot_mapped_reads,
        )
    else:
        logging.info("(Couldn't calculate total number of mapped reads)")

    return cnarr  # type: ignore[return-value]


def interval_coverages_bedgraph(
    regions_bed_fname: str,
    bedgraph_fname: str,
) -> CNA:
    """Calculate log2 coverages from bedGraph file at each interval.

    Parameters
    ----------
    regions_bed_fname : str
        Path to BED file defining regions to measure coverage.
    bedgraph_fname : str
        Path to tabix-indexed bedGraph file (.bed.gz).

    Returns
    -------
    CopyNumArray
        Coverage values for each region.
    """
    meta = {"sample_id": core.fbase(bedgraph_fname)}
    start_time = time.time()

    # Skip processing if the BED file is empty
    with open(regions_bed_fname) as bed_handle:
        for line in bed_handle:
            if line.strip():
                break
        else:
            logging.info(
                "Skip processing %s with empty regions file %s",
                os.path.basename(bedgraph_fname),
                regions_bed_fname,
            )
            return CNA.from_rows([], meta_dict=meta)  # type: ignore[return-value]

    # Calculate coverage from bedGraph
    logging.info("Processing bedGraph %s", os.path.basename(bedgraph_fname))
    table = bedgraph_to_basecount(regions_bed_fname, bedgraph_fname)

    # Fill in CNA required columns (same as interval_coverages_pileup)
    if "gene" in table:
        table["gene"] = table["gene"].fillna("-")
    else:
        table["gene"] = "-"

    # Calculate depth and log2 (same logic as interval_coverages_pileup)
    spans = table.end - table.start
    ok_idx = spans > 0
    table = table.assign(depth=0.0, log2=NULL_LOG2_COVERAGE)
    table.loc[ok_idx, "depth"] = table.loc[ok_idx, "basecount"] / spans[ok_idx]
    ok_idx = table["depth"] > 0
    table.loc[ok_idx, "log2"] = np.log2(table.loc[ok_idx, "depth"])

    # Remove basecount column before creating CNA
    table = table.drop("basecount", axis=1)
    cnarr = CNA(table, meta)

    # Log stats
    tot_time = time.time() - start_time
    logging.info(
        "Time: %.3f seconds (%d bins/sec)",
        tot_time,
        int(round(len(cnarr) / tot_time, 0)),
    )
    logging.info(
        "Summary: #bins=%d, mean_depth=%.4f, min_depth=%.4f, max_depth=%.4f",
        len(cnarr),
        cnarr["depth"].mean(),
        cnarr["depth"].min(),
        cnarr["depth"].max(),
    )

    return cnarr


def interval_coverages_count(
    bed_fname: str, bam_fname: str, min_mapq: int, procs: int = 1, fasta: None = None
) -> Iterator[list[int | tuple[str, int, int, str, float, float]]]:
    """Calculate log2 coverages in the BAM file at each interval."""
    try:
        import pysam
    except ImportError:
        raise ImportError(
            f"pysam is required for BAM read counting. {PYSAM_INSTALL_MSG}"
        ) from None
    regions = tabio.read_auto(bed_fname)
    # Skip regions on contigs absent from the BAM: pysam.fetch raises
    # "invalid contig" on them, same chrom-name-mismatch class as #620.
    bam_chroms = samutil.get_bam_chroms(bam_fname, fasta)
    present = [
        (chrom, subr) for chrom, subr in regions.by_chromosome() if chrom in bam_chroms
    ]
    if not present:
        raise ValueError(
            f"BED file {bed_fname!r} chromosome names don't match any in "
            f"BAM file {bam_fname!r}"
        )
    if procs == 1:
        bamfile = pysam.AlignmentFile(bam_fname, "rb", reference_filename=fasta)
        for chrom, subregions in present:
            logging.info(
                "Processing chromosome %s of %s", chrom, os.path.basename(bam_fname)
            )
            for count, row in _rdc_chunk(bamfile, subregions, min_mapq):
                yield [count, row]
    else:
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bam_fname, subr, min_mapq, fasta) for _c, subr in present)
            for chunk in pool.map(_rdc, args_iter):
                for count, row in chunk:
                    yield [count, row]


def _rdc(args):
    """Wrapper for parallel."""
    return list(_rdc_chunk(*args))


def _rdc_chunk(
    bamfile: pysam.AlignmentFile,
    regions: GenomicArray,
    min_mapq: int,
    fasta: None = None,
) -> Iterator[tuple[int, tuple[str, int, int, str, float, float]]]:
    if isinstance(bamfile, str):  # type: ignore[unreachable]
        try:  # type: ignore[unreachable]
            import pysam
        except ImportError:
            raise ImportError(
                f"pysam is required for BAM read counting. {PYSAM_INSTALL_MSG}"
            ) from None
        bamfile = pysam.AlignmentFile(bamfile, "rb", reference_filename=fasta)
    for chrom, start, end, gene in regions.coords(["gene"]):
        yield region_depth_count(bamfile, chrom, start, end, gene, min_mapq)


def region_depth_count(
    bamfile: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    gene: str,
    min_mapq: int,
) -> tuple[int, tuple[str, int, int, str, float, float]]:
    """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """

    def filter_read(read) -> bool:
        """True if the given read should be counted towards coverage.

        Excludes the same reads as ``COVERAGE_EXCLUDE_FLAGS`` (0x704) in the
        depth/bedcov path, plus low-MAPQ reads (#689).
        """
        return not (
            read.is_duplicate
            or read.is_secondary
            or read.is_unmapped
            or read.is_qcfail
            or read.mapq < min_mapq
        )

    count = 0
    bases = 0
    for read in bamfile.fetch(reference=chrom, start=start, end=end):
        if filter_read(read):
            count += 1
            # Only count the bases aligned to the region
            bases += sum(1 for p in read.positions if start <= p < end)  # type: ignore[attr-defined,misc]
    depth = bases / (end - start) if end > start else 0
    row = (
        chrom,
        start,
        end,
        gene,
        math.log(depth, 2) if depth else NULL_LOG2_COVERAGE,
        depth,
    )
    return count, row


def interval_coverages_pileup(
    bed_fname: str,
    bam_fname: str,
    min_mapq: int,
    procs: int = 1,
    fasta: str | None = None,
) -> pd.DataFrame:
    """Calculate log2 coverages in the BAM file at each interval."""
    logging.info("Processing reads in %s", os.path.basename(bam_fname))
    # Regions on contigs absent from the BAM header make samtools bedcov fail;
    # filter them out per call so a chrom-name mismatch on some (but not all)
    # regions doesn't abort the run -- and so the result is independent of how
    # the BED is split across processes (#620).
    bam_chroms = samutil.get_bam_chroms(bam_fname, fasta)
    if procs == 1:
        table = bedcov(bed_fname, bam_fname, min_mapq, fasta, bam_chroms=bam_chroms)
    else:
        chunks = []
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = (
                (bed_chunk, bam_fname, min_mapq, fasta, bam_chroms)
                for bed_chunk in to_chunks(bed_fname)
            )
            for bed_chunk_fname, table in pool.map(_bedcov, args_iter):
                if len(table):
                    chunks.append(table)
                rm(bed_chunk_fname)
        if not chunks:
            raise ValueError(
                f"BED file {bed_fname!r} chromosome names don't match any in "
                f"BAM file {bam_fname!r}"
            )
        table = pd.concat(chunks, ignore_index=True)
    # Fill in CNA required columns
    if "gene" in table:
        table["gene"] = table["gene"].fillna("-")
    else:
        table["gene"] = "-"
    # User-supplied bins might be zero-width or reversed -- skip those
    spans = table.end - table.start
    ok_idx = spans > 0
    table = table.assign(depth=0.0, log2=NULL_LOG2_COVERAGE)
    table.loc[ok_idx, "depth"] = table.loc[ok_idx, "basecount"] / spans[ok_idx]
    ok_idx = table["depth"] > 0
    table.loc[ok_idx, "log2"] = np.log2(table.loc[ok_idx, "depth"])
    return table


def _bedcov(args):
    """Wrapper for parallel: filter each chunk to BAM contigs, tolerate empties.

    Per-chunk emptiness is expected (a chunk may land entirely on contigs
    absent from the BAM); the caller raises if *all* chunks come back empty.
    """
    bed_fname, bam_fname, min_mapq, fasta, bam_chroms = args
    table = bedcov(
        bed_fname, bam_fname, min_mapq, fasta, allow_empty=True, bam_chroms=bam_chroms
    )
    return bed_fname, table


def _filter_bed_to_chroms(bed_fname: str, keep_chroms: set[str]) -> tuple[str, int]:
    """Copy `bed_fname` keeping only lines on chromosomes in `keep_chroms`.

    Returns the path to a new temporary BED file and the number of lines kept.
    The caller is responsible for removing the temporary file.
    """
    fd, tmp_path = tempfile.mkstemp(suffix=".bed", prefix="cnvkit-filt.")
    opener = gzip.open if bed_fname.endswith(".gz") else open
    n_kept = 0
    # fdopen first so the temp fd is always closed even if opening the input
    # raises (e.g. a malformed .gz).
    with os.fdopen(fd, "w") as out, opener(bed_fname, "rt") as infile:
        for line in infile:
            if not line.strip() or line[0] == "#":
                continue
            if line.split("\t", 1)[0] in keep_chroms:
                out.write(line)
                n_kept += 1
    return tmp_path, n_kept


def bedcov(
    bed_fname: str,
    bam_fname: str,
    min_mapq: int,
    fasta: str | None = None,
    max_depth: int | None = None,
    allow_empty: bool = False,
    bam_chroms: set[str] | None = None,
) -> pd.DataFrame:
    """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.

    i.e. mean pileup depth across each region.

    Parameters
    ----------
    bed_fname : str
        Path to BED file defining regions to measure.
    bam_fname : str
        Path to the indexed BAM/CRAM file.
    min_mapq : int
        Minimum read mapping quality to count.
    fasta : str, optional
        Reference FASTA, required for CRAM input.
    max_depth : int, optional
        Cap pileup depth at this value (samtools -d).
    allow_empty : bool
        If True, return an empty DataFrame instead of raising when no regions
        remain to measure (e.g. a parallel chunk landing entirely on contigs
        absent from the BAM). Caller validates global chrom-name concordance.
    bam_chroms : set of str, optional
        Reference names present in the BAM header. When given, BED regions on
        any other contig are dropped before calling samtools, which otherwise
        aborts on regions whose contig is absent from the BAM (#620).
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            f"pysam is required for BAM coverage calculation. {PYSAM_INSTALL_MSG}"
        ) from None
    # Drop regions on contigs absent from the BAM so samtools doesn't abort
    filtered_fname = None
    if bam_chroms is not None:
        filtered_fname, n_kept = _filter_bed_to_chroms(bed_fname, bam_chroms)
        if not n_kept:
            rm(filtered_fname)
            if allow_empty:
                return pd.DataFrame()
            raise ValueError(
                f"BED file {bed_fname!r} chromosome names don't match any in "
                f"BAM file {bam_fname!r}"
            )
    try:
        # Count bases in each region; exclude low-MAPQ and flagged reads.
        # `-G` *adds* to bedcov's default exclude set, so this guarantees
        # unmapped/secondary/QC-fail/duplicate reads are skipped regardless of
        # the samtools default (#689).
        cmd = ["-G", hex(COVERAGE_EXCLUDE_FLAGS)]
        if max_depth is not None:
            cmd.extend(["-d", str(max_depth)])
        if min_mapq and min_mapq > 0:
            cmd.extend(["-Q", str(min_mapq)])
        if fasta:
            cmd.extend(["--reference", fasta])
        cmd.extend([filtered_fname or bed_fname, bam_fname])
        try:
            raw = pysam.bedcov(*cmd, split_lines=False)  # type: ignore[attr-defined]
        except pysam.SamtoolsError as exc:
            raise ValueError(
                f"Failed processing {bam_fname!r} coverages in {bed_fname!r} "
                f"regions. PySAM error: {exc}"
            ) from exc
    finally:
        if filtered_fname is not None:
            rm(filtered_fname)
    if not raw:
        if allow_empty:
            return pd.DataFrame()
        raise ValueError(
            f"BED file {bed_fname!r} chromosome names don't match any in "
            f"BAM file {bam_fname!r}"
        )
    columns = detect_bedcov_columns(raw)  # type: ignore[arg-type]
    usecols = [c for c in columns if c != "extra"]
    table = pd.read_csv(StringIO(raw), sep="\t", names=columns, usecols=usecols)  # type: ignore[arg-type]
    return table


def detect_bedcov_columns(text: str) -> list[str]:
    """Determine which 'bedcov' output columns to keep.

    bedcov outputs the input BED columns plus an appended numeric column
    (basecount). With some options (e.g. -d), an additional trailing numeric
    column may be appended; then basecount is the second-to-last column.
    """
    firstline = text[: text.index("\n")]
    fields = firstline.split("\t")
    tabcount = len(fields) - 1

    if tabcount < 3:
        raise RuntimeError(f"Bad line from bedcov:\n{firstline!r}")

    def _is_int(s: str) -> bool:
        try:
            int(s)
            return True
        except ValueError:
            return False

    has_extra = len(fields) >= 2 and _is_int(fields[-1]) and _is_int(fields[-2])

    # BED3
    if tabcount == 3:
        return (
            ["chromosome", "start", "end", "basecount", "extra"]
            if has_extra
            else ["chromosome", "start", "end", "basecount"]
        )

    # BED4
    if tabcount == 4:
        return (
            ["chromosome", "start", "end", "gene", "basecount", "extra"]
            if has_extra
            else ["chromosome", "start", "end", "gene", "basecount"]
        )

    # BED4+ with extra columns after gene
    # Input BED has arbitrary columns after 'gene' -- ignore them
    n_numeric = 2 if has_extra else 1
    n_total = len(fields)  # total columns in output
    n_fillers = n_total - 4 - n_numeric

    if n_fillers < 0:
        raise RuntimeError(f"Unexpected bedcov output:\n{firstline!r}")

    fillers = [f"_{i}" for i in range(1, n_fillers + 1)]
    cols = ["chromosome", "start", "end", "gene", *fillers, "basecount"]
    if has_extra:
        cols.append("extra")
    return cols
