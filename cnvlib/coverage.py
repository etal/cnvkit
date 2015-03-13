"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division
import math
import os.path
import time
from itertools import groupby

from Bio._py3k import map, zip
import pysam

from . import ngfrills
from .ngfrills import echo

from .params import NULL_LOG2_COVERAGE, READ_LEN


def interval_coverages(bed_fname, bam_fname, by_count):
    """Calculate log2 coverages in the BAM file at each interval."""
    start_time = time.time()

    # Skip processing if the BED file is empty
    with open(bed_fname) as bed_handle:
        for line in bed_handle:
            if line.strip():
                break
        else:
            echo("Skip processing", os.path.basename(bam_fname),
                "with empty regions file", bed_fname)
            return []

    # Calculate average read depth in each bin
    ic_func = (interval_coverages_count if by_count
               else interval_coverages_pileup)
    results = ic_func(bed_fname, bam_fname)
    read_counts, cna_rows = zip(*list(results))

    # Log some stats
    tot_time = time.time() - start_time
    tot_reads = sum(read_counts)
    echo("Time: %.3f seconds (%d reads/sec, %s bins/sec)"
         % (tot_time,
            int(round(tot_reads / tot_time, 0)),
            int(round(len(read_counts) / tot_time, 0))))
    echo("Summary:",
         "#bins=%d," % len(read_counts),
         "#reads=%d," % tot_reads,
         "mean=%.4f," % (tot_reads / len(read_counts)),
         "min=%s," % min(read_counts),
         "max=%s" % max(read_counts))
    tot_mapped_reads = bam_total_reads(bam_fname)
    if tot_mapped_reads:
        echo("On-target percentage: %.3f (of %d mapped)"
            % (100. * tot_reads / tot_mapped_reads, tot_mapped_reads))
    else:
        echo("(Couldn't calculate total number of mapped reads)")

    return list(cna_rows)


def interval_coverages_count(bed_fname, bam_fname):
    """Calculate log2 coverages in the BAM file at each interval."""
    bamfile = pysam.Samfile(bam_fname, 'rb')
    # Parse the BED lines and group them by chromosome
    # (efficient if records are already sorted by chromosome)
    for chrom, rows_iter in groupby(ngfrills.parse_regions(bed_fname),
                                    key=lambda r: r[0]):
        # Thunk and reshape this chromosome's intervals
        echo("Processing chromosome", chrom, "of", os.path.basename(bam_fname))
        _chroms, starts, ends, names = zip(*rows_iter)
        counts_depths = [region_depth_count(bamfile, chrom, s, e)
                         for s, e in zip(starts, ends)]
        for start, end, name, (count, depth) in zip(starts, ends, names,
                                                    counts_depths):
            yield [count,
                   (chrom, start, end, name,
                    math.log(depth, 2) if depth else NULL_LOG2_COVERAGE)]


def region_depth_count(bamfile, chrom, start, end):
    """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """
    # Count the number of read midpoints in the interval
    count = sum((start <= (read.pos + .5*read.rlen) <= end and filter_read(read))
                for read in bamfile.fetch(reference=chrom,
                                          start=start, end=end))
    # Scale read counts to region length
    # Depth := #Bases / Span
    depth = (READ_LEN * count / (end - start)
             if end > start else 0)
    return count, depth


def interval_coverages_pileup(bed_fname, bam_fname):
    """Calculate log2 coverages in the BAM file at each interval."""
    echo("Processing reads in", os.path.basename(bam_fname))
    for chrom, start, end, name, count, depth in bedcov(bed_fname, bam_fname):
        yield [count,
               (chrom, start, end, name,
                math.log(depth, 2) if depth else NULL_LOG2_COVERAGE)]


def bedcov(bed_fname, bam_fname):
    """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.

    i.e. mean pileup depth across each region.
    """
    # Count bases in each region; exclude 0-MAPQ reads
    try:
        lines = pysam.bedcov(bed_fname, bam_fname, '-Q', '1')
    except pysam.SamtoolsError as exc:
        raise ValueError("Failed processing %r coverages in %r regions. PySAM error: %s"
                         % (bam_fname, bed_fname, exc))
    if not lines:
        raise ValueError("BED file %r sequence IDs don't match any in BAM file %r"
                         % (bed_fname, bam_fname))
    # Return an iterable...
    for line in lines:
        try:
            chrom, start_s, end_s, name, basecount_s = line.split('\t')
        except:
            raise RuntimeError("Bad line from bedcov:\n" + line)
        start, end, basecount = map(int, (start_s, end_s, basecount_s.strip()))
        span = end - start
        if span > 0:
            # Algebra from above
            count = basecount / READ_LEN
            mean_depth = basecount / span
        else:
            # User-supplied bins might be oddly constructed
            count = mean_depth = 0
        yield chrom, start, end, name, count, mean_depth


def filter_column(col):
    """Count the number of filtered reads in a pileup column."""
    return sum(filter_read(read.alignment) for read in col.pileups)


def filter_read(read):
    """True if the given read should be counted towards coverage."""
    return not (read.is_duplicate
                or read.is_secondary
                or read.is_unmapped
                or read.is_qcfail
                or read.mapq == 0
               )


def bam_total_reads(bam_fname):
    """Count the total number of mapped reads in a BAM file.

    Uses the BAM index to do this quickly.
    """
    lines = pysam.idxstats(bam_fname)
    tot_mapped_reads = 0
    for line in lines:
        _seqname, _seqlen, nmapped, _nunmapped = line.split()
        tot_mapped_reads += int(nmapped)
    return tot_mapped_reads

