"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division, print_function
from builtins import zip
from past.builtins import basestring

import logging
import math
import os.path
import time
from concurrent import futures

import numpy as np
import pandas as pd
import pysam
from Bio._py3k import StringIO
from skgenome import tabio

from . import core, samutil
from .cnary import CopyNumArray as CNA
from .parallel import rm, to_chunks
from .params import NULL_LOG2_COVERAGE


def do_coverage(bed_fname, bam_fname, by_count=False, min_mapq=0, processes=1):
    """Calculate coverage in the given regions from BAM read depths."""
    if not samutil.ensure_bam_sorted(bam_fname):
        raise RuntimeError("BAM file %s must be sorted by coordinates"
                            % bam_fname)
    samutil.ensure_bam_index(bam_fname)
    # ENH: count importers.TOO_MANY_NO_COVERAGE & warn
    cnarr = interval_coverages(bed_fname, bam_fname, by_count, min_mapq,
                               processes)
    return cnarr


def interval_coverages(bed_fname, bam_fname, by_count, min_mapq, processes):
    """Calculate log2 coverages in the BAM file at each interval."""
    meta = {'sample_id': core.fbase(bam_fname)}
    start_time = time.time()

    # Skip processing if the BED file is empty
    with open(bed_fname) as bed_handle:
        for line in bed_handle:
            if line.strip():
                break
        else:
            logging.info("Skip processing %s with empty regions file %s",
                         os.path.basename(bam_fname), bed_fname)
            return CNA.from_rows([], meta_dict=meta)

    # Calculate average read depth in each bin
    if by_count:
        results = interval_coverages_count(bed_fname, bam_fname, min_mapq,
                                           processes)
        read_counts, cna_rows = zip(*results)
        read_counts = pd.Series(read_counts)
        cnarr = CNA.from_rows(list(cna_rows),
                              columns=CNA._required_columns + ('depth',),
                              meta_dict=meta)
    else:
        table = interval_coverages_pileup(bed_fname, bam_fname, min_mapq,
                                          processes)
        read_len = samutil.get_read_length(bam_fname)
        read_counts = table['basecount'] / read_len
        table = table.drop('basecount', axis=1)
        cnarr = CNA(table, meta)

    # Log some stats
    tot_time = time.time() - start_time
    tot_reads = read_counts.sum()
    logging.info("Time: %.3f seconds (%d reads/sec, %s bins/sec)",
                 tot_time,
                 int(round(tot_reads / tot_time, 0)),
                 int(round(len(read_counts) / tot_time, 0)))
    logging.info("Summary: #bins=%d, #reads=%d, "
                 "mean=%.4f, min=%s, max=%s ",
                 len(read_counts),
                 tot_reads,
                 (tot_reads / len(read_counts)),
                 read_counts.min(),
                 read_counts.max())
    tot_mapped_reads = samutil.bam_total_reads(bam_fname)
    if tot_mapped_reads:
        logging.info("Percent reads in regions: %.3f (of %d mapped)",
                     100. * tot_reads / tot_mapped_reads,
                     tot_mapped_reads)
    else:
        logging.info("(Couldn't calculate total number of mapped reads)")

    return cnarr


def interval_coverages_count(bed_fname, bam_fname, min_mapq, procs=1):
    """Calculate log2 coverages in the BAM file at each interval."""
    regions = tabio.read_auto(bed_fname)
    if procs == 1:
        bamfile = pysam.Samfile(bam_fname, 'rb')
        for chrom, subregions in regions.by_chromosome():
            logging.info("Processing chromosome %s of %s",
                         chrom, os.path.basename(bam_fname))
            for count, row in _rdc_chunk(bamfile, subregions, min_mapq):
                yield [count, row]
    else:
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bam_fname, subr, min_mapq)
                         for _c, subr in regions.by_chromosome())
            for chunk in pool.map(_rdc, args_iter):
                for count, row in chunk:
                    yield [count, row]


def _rdc(args):
    """Wrapper for parallel."""
    return list(_rdc_chunk(*args))


def _rdc_chunk(bamfile, regions, min_mapq):
    if isinstance(bamfile, basestring):
        bamfile = pysam.Samfile(bamfile, 'rb')
    for chrom, start, end, gene in regions.coords(["gene"]):
        yield region_depth_count(bamfile, chrom, start, end, gene, min_mapq)


def region_depth_count(bamfile, chrom, start, end, gene, min_mapq):
    """Calculate depth of a region via pysam count.

    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.

    Coordinates are 0-based, per pysam.
    """
    def filter_read(read):
        """True if the given read should be counted towards coverage."""
        return not (read.is_duplicate
                    or read.is_secondary
                    or read.is_unmapped
                    or read.is_qcfail
                    or read.mapq < min_mapq)

    count = 0
    bases = 0
    for read in bamfile.fetch(reference=chrom, start=start, end=end):
        if filter_read(read):
            count += 1
            # Only count the bases aligned to the region
            rlen = read.query_length
            if read.pos < start:
                rlen -= start - read.pos
            if read.pos + read.query_length > end:
                rlen -= read.pos + read.query_length - end
            bases += rlen
    depth = bases / (end - start) if end > start else 0
    row = (chrom, start, end, gene,
           math.log(depth, 2) if depth else NULL_LOG2_COVERAGE,
           depth)
    return count, row


def interval_coverages_pileup(bed_fname, bam_fname, min_mapq, procs=1):
    """Calculate log2 coverages in the BAM file at each interval."""
    logging.info("Processing reads in %s", os.path.basename(bam_fname))
    if procs == 1:
        table = bedcov(bed_fname, bam_fname, min_mapq)
    else:
        chunks = []
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bed_chunk, bam_fname, min_mapq)
                         for bed_chunk in to_chunks(bed_fname))
            for bed_chunk_fname, table in pool.map(_bedcov, args_iter):
                chunks.append(table)
                rm(bed_chunk_fname)
        table = pd.concat(chunks, ignore_index=True)
    # Fill in CNA required columns
    if 'gene' in table:
        table['gene'] = table['gene'].fillna('-')
    else:
        table['gene'] = '-'
    # User-supplied bins might be zero-width or reversed -- skip those
    spans = table.end - table.start
    ok_idx = (spans > 0)
    table = table.assign(depth=0, log2=NULL_LOG2_COVERAGE)
    table.loc[ok_idx, 'depth'] = (table.loc[ok_idx, 'basecount']
                                  / spans[ok_idx])
    ok_idx = (table['depth'] > 0)
    table.loc[ok_idx, 'log2'] = np.log2(table.loc[ok_idx, 'depth'])
    return table


def _bedcov(args):
    """Wrapper for parallel."""
    bed_fname = args[0]
    table = bedcov(*args)
    return bed_fname, table


def bedcov(bed_fname, bam_fname, min_mapq):
    """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.

    i.e. mean pileup depth across each region.
    """
    # Count bases in each region; exclude low-MAPQ reads
    cmd = [bed_fname, bam_fname]
    if min_mapq and min_mapq > 0:
        cmd.extend(['-Q', bytes(min_mapq)])
    try:
        raw = pysam.bedcov(*cmd, split_lines=False)
    except pysam.SamtoolsError as exc:
        raise ValueError("Failed processing %r coverages in %r regions. "
                         "PySAM error: %s" % (bam_fname, bed_fname, exc))
    if not raw:
        raise ValueError("BED file %r chromosome names don't match any in "
                         "BAM file %r" % (bed_fname, bam_fname))
    columns = detect_bedcov_columns(raw)
    table = pd.read_table(StringIO(raw), names=columns, usecols=columns)
    return table


def detect_bedcov_columns(text):
    """Determine which 'bedcov' output columns to keep.

    Format is the input BED plus a final appended column with the count of
    basepairs mapped within each row's region. The input BED might have 3
    columns (regions without names), 4 (named regions), or more (arbitrary
    columns after 'gene').
    """
    firstline = text[:text.index('\n')]
    tabcount = firstline.count('\t')
    if tabcount < 3:
        raise RuntimeError("Bad line from bedcov:\n%r" % firstline)
    if tabcount == 3:
        return ['chromosome', 'start', 'end', 'basecount']
    if tabcount == 4:
        return ['chromosome', 'start', 'end', 'gene', 'basecount']
    # Input BED has arbitrary columns after 'gene' -- ignore them
    fillers = ["_%d" % i for i in range(1, tabcount - 3)]
    return ['chromosome', 'start', 'end', 'gene'] + fillers + ['basecount']
