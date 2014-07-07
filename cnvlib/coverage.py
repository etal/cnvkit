"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division
import math
import os.path
import time
from itertools import groupby

from Bio._py3k import map, range, zip
import pysam

from . import ngfrills
from .ngfrills import echo

from .params import NULL_LOG2_COVERAGE, READ_LEN


def interval_coverages(bed_fname, bam_fname, region_depth_func):
    """Calculate log2 coverages in the BAM file at each interval."""
    start_time = time.time()

    bamfile = pysam.Samfile(bam_fname, 'rb')
    # Parse the BED lines and group them by chromosome
    # (efficient if records are already sorted by chromosome)
    all_counts = []
    for chrom, rows_iter in groupby(ngfrills.parse_regions(bed_fname),
                                    key=lambda r: r[0]):
        # Thunk and reshape this chromosome's intervals
        echo("Processing chromosome", chrom, "of", os.path.basename(bam_fname))
        _chroms, starts, ends, names = zip(*rows_iter)
        counts_depths = [region_depth_func(bamfile, chrom, s, e)
                         for s, e in zip(starts, ends)]
        for start, end, name, (count, depth) in zip(starts, ends, names,
                                                    counts_depths):
            all_counts.append(count)
            yield (chrom, start, end, name,
                   math.log(depth, 2) if depth else NULL_LOG2_COVERAGE)
    # Log some stats
    tot_time = time.time() - start_time
    tot_reads = sum(all_counts)
    echo("Time: %.3f seconds (%d reads/sec, %s bins/sec)"
         % (tot_time,
            int(round(tot_reads / tot_time, 0)),
            int(round(len(all_counts) / tot_time, 0))))
    echo("Summary of counts:",
         "\n\tbins=%d, total reads=%d" % (len(all_counts), tot_reads),
         "\n\tmean=%.4f min=%s max=%s" % (tot_reads / len(all_counts),
                                          min(all_counts), max(all_counts)))
    tot_mapped_reads = bam_total_reads(bam_fname)
    if tot_mapped_reads:
        echo("On-target percentage: %.3f (of %d mapped)"
            % (100. * tot_reads / tot_mapped_reads, tot_mapped_reads))
    else:
        echo("(Couldn't calculate total number of mapped reads)")


def region_depth_count(bamfile, chrom, start, end):
    """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """
    # ENH: shrink/stretch region by average read length?
    # Count the number of read start positions in the interval (like Picard)
    count = sum((start <= read.pos <= end and filter_read(read))
                for read in bamfile.fetch(reference=chrom,
                                          start=start, end=end))
    # Scale read counts to region length
    # Depth := #Bases / Span
    depth = READ_LEN * count / (end - start)
    return count, depth


def region_depth_pileup(bamfile, chrom, start, end):
    """Calculate depth of a region via pysam pileup.

    i.e. average pileup depth across a region.

    To reduce edge effects, depth is counted within a narrowed region, shifting
    the start and end points inward by a margin equal to the expected read
    length. (Or, if the region is too narrow, take the middle two quartiles.)

    Coordinates are 0-based, per pysam.
    """
    # Narrow the region to reduce edge effects
    # quartile = .25 * (end - start)
    # start = min(start + READ_LEN, start + quartile)
    # end = max(end - READ_LEN, end - quartile)
    depths = [filter_column(col) for col in bamfile.pileup(chrom, start, end)]
    mean_depth = sum(depths) / len(depths) if depths else 0
    count = mean_depth * (end - start) / READ_LEN  # algebra from above
    return count, mean_depth  # Mean


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

