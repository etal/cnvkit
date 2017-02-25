"""BAM utilities."""
from __future__ import absolute_import, division, print_function
#  from builtins import str
from past.builtins import basestring

import logging
import os
from itertools import islice

import numpy as np
import pandas as pd
import pysam
from Bio._py3k import StringIO


def idxstats(bam_fname, drop_unmapped=False):
    """Get chromosome names, lengths, and number of mapped/unmapped reads.

    Use the BAM index (.bai) to get the number of reads and size of each
    chromosome. Contigs with no mapped reads are skipped.
    """
    handle = StringIO(pysam.idxstats(bam_fname, split_lines=False))
    table = pd.read_table(handle, header=None,
                          names=['chromosome', 'length', 'mapped', 'unmapped'])
    if drop_unmapped:
        table = table[table.mapped != 0].drop('unmapped', axis=1)
    return table


def bam_total_reads(bam_fname):
    """Count the total number of mapped reads in a BAM file.

    Uses the BAM index to do this quickly.
    """
    table = idxstats(bam_fname, drop_unmapped=True)
    return table.mapped.sum()


def ensure_bam_index(bam_fname):
    """Ensure a BAM file is indexed, to enable fast traversal & lookup.

    For MySample.bam, samtools will look for an index in these files, in order:

    - MySample.bam.bai
    - MySample.bai
    """
    if os.path.isfile(bam_fname + '.bai'):
        # MySample.bam.bai
        bai_fname = bam_fname + '.bai'
    else:
        # MySample.bai
        bai_fname = bam_fname[:-1] + 'i'
    if not is_newer_than(bai_fname, bam_fname):
        logging.info("Indexing BAM file %s", bam_fname)
        pysam.index(bam_fname)
        bai_fname = bam_fname + '.bai'
    assert os.path.isfile(bai_fname), \
            "Failed to generate index " + bai_fname
    return bai_fname


def ensure_bam_sorted(bam_fname, by_name=False, span=50):
    """Test if the reads in a BAM file are sorted as expected.

    by_name=True: reads are expected to be sorted by query name. Consecutive
    read IDs are in alphabetical order, and read pairs appear together.

    by_name=False: reads are sorted by position. Consecutive reads have
    increasing position.
    """
    if by_name:
        # Compare read IDs
        def out_of_order(read, prev):
            return not (prev is None or
                        prev.qname <= read.qname)
    else:
        # Compare read locations
        def out_of_order(read, prev):
            return not (prev is None or
                        read.tid != prev.tid or
                        prev.pos <= read.pos)

    # ENH - repeat at 50%, ~99% through the BAM
    bam = pysam.Samfile(bam_fname, 'rb')
    last_read = None
    for read in islice(bam, span):
        if out_of_order(read, last_read):
            return False
        last_read = read
    bam.close()
    return True


def is_newer_than(target_fname, orig_fname):
    """Compare file modification times."""
    if not os.path.isfile(target_fname):
        return False
    return (os.stat(target_fname).st_mtime >= os.stat(orig_fname).st_mtime)


def get_read_length(bam, span=1000):
    """Get (median) read length from first few reads in a BAM file.

    Illumina reads all have the same length; other sequencers might not.

    Parameters
    ----------
    bam : str or pysam.Samfile
        Filename or pysam-opened BAM file.
    n : int
        Number of reads used to calculate median read length.
    """
    was_open = False
    if isinstance(bam, basestring):
        bam = pysam.Samfile(bam, 'rb')
    else:
        was_open = True
    lengths = [read.query_length for read in islice(bam, span)
               if read.query_length > 0]
    if was_open:
        bam.seek(0)
    else:
        bam.close()
    return np.median(lengths)
