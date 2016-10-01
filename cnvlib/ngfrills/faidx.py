"""NGS utilities: Indexed FASTA I/O."""
from __future__ import absolute_import, division, print_function
from builtins import str

import logging
from pyfaidx import Fasta


def fasta_extract_regions(fa_fname, intervals):
    """Extract an iterable of regions from an indexed FASTA file.

    Input: FASTA file name; iterable of (seq_id, start, end) (1-based)
    Output: iterable of string sequences.
    """
    with Fasta(fa_fname, as_raw=True) as fa_file:
        for chrom, subarr in intervals.by_chromosome():
            logging.info("Extracting sequences from chromosome %s", chrom)
            for _chrom, start, end in subarr.coords():
                yield fa_file[_chrom][start.item():end.item()]


def _fasta_extract_regions_safe(fa_fname, intervals):
    """Simpler, slower version of fasta_extract_regions, for testing it."""
    from Bio import SeqIO
    idx = SeqIO.index(fa_fname, 'fasta')
    for chrom, subarr in intervals.by_chromosome():
        logging.info("Extracting sequences from chromosome %s", chrom)
        seq = str(idx[chrom].seq)
        for _chrom, start, end in subarr.coords():
            subseq = seq[start:end]
            if len(subseq) != end - start:
                logging.info("Short subsequence %s:%d-%d; read %d, wanted %d",
                             _chrom, start, end, len(subseq), end - start)
            yield subseq
