"""Supporting functions for the 'antitarget' command."""
from __future__ import absolute_import, division, print_function
from builtins import str, zip
from past.builtins import basestring

import atexit
import gzip
import logging
import math
import os
import os.path
import tempfile
import time
from concurrent import futures
from io import StringIO

import numpy as np
import pandas as pd
import pysam

from . import tabio
from .cnary import CopyNumArray as CNA
from .core import fbase
from .params import NULL_LOG2_COVERAGE, READ_LEN
from .samutil import bam_total_reads


def interval_coverages(bed_fname, bam_fname, by_count, min_mapq, processes):
    """Calculate log2 coverages in the BAM file at each interval."""
    start_time = time.time()

    # Skip processing if the BED file is empty
    with open(bed_fname) as bed_handle:
        for line in bed_handle:
            if line.strip():
                break
        else:
            logging.info("Skip processing %s with empty regions file %s",
                         os.path.basename(bam_fname), bed_fname)
            return CNA.from_rows([], meta_dict={'sample_id': fbase(bam_fname)})

    # Calculate average read depth in each bin
    ic_func = (interval_coverages_count if by_count
               else interval_coverages_pileup)
    results = ic_func(bed_fname, bam_fname, min_mapq, processes)
    read_counts, cna_rows = list(zip(*list(results)))

    # Log some stats
    tot_time = time.time() - start_time
    tot_reads = sum(read_counts)
    logging.info("Time: %.3f seconds (%d reads/sec, %s bins/sec)",
                 tot_time,
                 int(round(tot_reads / tot_time, 0)),
                 int(round(len(read_counts) / tot_time, 0)))
    logging.info("Summary: #bins=%d, #reads=%d, "
                 "mean=%.4f, min=%s, max=%s ",
                 len(read_counts),
                 tot_reads,
                 (tot_reads / len(read_counts)),
                 min(read_counts),
                 max(read_counts))
    tot_mapped_reads = bam_total_reads(bam_fname)
    if tot_mapped_reads:
        logging.info("Percent reads in regions: %.3f (of %d mapped)",
                     100. * tot_reads / tot_mapped_reads,
                     tot_mapped_reads)
    else:
        logging.info("(Couldn't calculate total number of mapped reads)")

    return CNA.from_rows(list(cna_rows),
                         columns=CNA._required_columns + ('depth',),
                         meta_dict={'sample_id': fbase(bam_fname)})


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

    # Count the number of read midpoints in the interval
    count = sum((start <= (read.pos + .5*read.rlen) <= end and filter_read(read))
                for read in bamfile.fetch(reference=chrom,
                                          start=start, end=end))
    # Scale read counts to region length
    # Depth := #Bases / Span
    depth = (READ_LEN * count / (end - start)
             if end > start else 0)
    row = (chrom, start, end, gene,
           math.log(depth, 2) if depth else NULL_LOG2_COVERAGE, depth)
    return count, row


def rm(path):
    """Safely remove a file."""
    try:
        os.unlink(path)
    except OSError:
        pass


def to_chunks(bed_fname, chunk_size=5000):
    """Split the bed-file into chunks for parallelization"""
    k, chunk = 0, 0
    fd, name = tempfile.mkstemp(suffix=".bed", prefix="tmp.%s." % chunk)
    fh = os.fdopen(fd, "w")
    atexit.register(rm, name)
    for line in (gzip.open if bed_fname.endswith(".gz") else open)(bed_fname):
        if line[0] == "#":
            continue
        k += 1
        fh.write(line)
        if k % chunk_size == 0:
            fh.close()
            yield name
            chunk += 1
            fd, name = tempfile.mkstemp(suffix=".bed", prefix="tmp.%s." % chunk)
            fh = os.fdopen(fd, "w")
    fh.close()
    if k % chunk_size:
        fh.close()
        yield name


def interval_coverages_pileup(bed_fname, bam_fname, min_mapq, procs=1):
    """Calculate log2 coverages in the BAM file at each interval."""
    logging.info("Processing reads in %s", os.path.basename(bam_fname))
    if procs == 1:
        for count, row in bedcov(bed_fname, bam_fname, min_mapq):
            yield [count, row]
    else:
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bed_chunk, bam_fname, min_mapq)
                         for bed_chunk in to_chunks(bed_fname))
            for bed_chunk_fname, biter in pool.map(_bedcov, args_iter):
                for count, row in biter:
                    yield [count, row]
                rm(bed_chunk_fname)


def _bedcov(args):
    """Wrapper for parallel."""
    bed_fname = args[0]
    return bed_fname, list(bedcov(*args))


def bedcov(bed_fname, bam_fname, min_mapq):
    """Calculate depth of all regions in a BED file via samtools (pysam) bedcov.

    i.e. mean pileup depth across each region.
    """
    # Count bases in each region; exclude low-MAPQ reads
    cmd = [bed_fname, bam_fname]
    if min_mapq and min_mapq > 0:
        cmd.extend(['-Q', bytes(min_mapq)])
    try:
        raw = pysam.bedcov(*cmd)
    except pysam.SamtoolsError as exc:
        raise ValueError("Failed processing %r coverages in %r regions. "
                         "PySAM error: %s" % (bam_fname, bed_fname, exc))
    if not raw:
        raise ValueError("BED file %r chromosome names don't match any in "
                         "BAM file %r" % (bed_fname, bam_fname))
    columns = detect_bedcov_columns(raw)
    table = pd.read_table(StringIO(str(raw)), names=columns, usecols=columns)
    if 'gene' in table:
        table['gene'] = table['gene'].fillna('-')
    else:
        table['gene'] = '-'
    spans = table.end - table.start
    # NB: User-supplied bins might be zero-width or reversed -- skip those
    ok_idx = (spans > 0)
    table = table.assign(count=0, mean_depth=0, log2=NULL_LOG2_COVERAGE)
    table.loc[ok_idx, 'mean_depth'] = (table.loc[ok_idx, 'basecount']
                                       / spans[ok_idx])
    table.loc[ok_idx, 'log2'] = np.log2(table.loc[ok_idx, 'mean_depth'])
    table.loc[ok_idx, 'count'] = table.loc[ok_idx, 'basecount'] / READ_LEN
    # Return an iterable...
    # ENH: return dframe here and above; pd.concat
    #   skip emitting count; calculate as: (depth * span)/READ_LEN
    for row in table.itertuples(index=False):
        yield row.count, (row.chromosome,
                          row.start,
                          row.end,
                          row.gene,
                          row.log2,
                          row.mean_depth)


def detect_bedcov_columns(text):
    # NB: gene is optional; may be trailing cols
    firstline = text[:text.index('\n')]
    tabcount = firstline.count('\t')
    if tabcount >= 4:
        return ['chromosome', 'start', 'end', 'gene', 'basecount']
    elif tabcount == 3:
        return ['chromosome', 'start', 'end', 'basecount']
    raise RuntimeError("Bad line from bedcov:\n%r" % firstline)
