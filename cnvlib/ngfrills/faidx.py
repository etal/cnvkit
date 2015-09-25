"""NGS utilities: Indexed FASTA I/O."""
from __future__ import absolute_import, division, print_function

import math
import os
from itertools import groupby

import pysam
from Bio._py3k import map

from .shared import echo, is_newer_than


def ensure_fasta_index(fasta_fname):
    """Ensure a FASTA file is indexed for samtools, to enable fast lookup."""
    fai_fname = fasta_fname + '.fai'
    if not is_newer_than(fai_fname, fasta_fname):
        echo("Indexing FASTA file", fasta_fname)
        pysam.faidx(fasta_fname)
    assert os.path.isfile(fai_fname), "Failed to generate index " + fai_fname
    return fai_fname


def read_fasta_index(fasta_fname):
    """Load a FASTA file's index.

    Returns a dict of:
        {seq_id: (length, offset, chars_per_line, bytes_per_line), ...}

    The index file contains, in one row per sequence, tab-separated columns:

        - sequence identifier
        - length
        - offset of the first sequence character in the file
        - number of characters per line
        - number of bytes per line (including the end-of-line character)

    With this information, we can easily compute the byte offset of the i-th
    character of a sequence in a file by looking at its index record. We skip to
    this byte offset in the file and from there, we can read the necessary
    sequence characters.

    See: http://trac.seqan.de/wiki/Tutorial/IndexedFastaIO:
    """
    index = {}
    fai_fname = ensure_fasta_index(fasta_fname)
    with open(fai_fname) as faifile:
        for line in faifile:
            fields = line.rstrip().split('\t')
            seq_id = fields[0]
            assert seq_id not in index, "Duplicate ID: " + seq_id
            index[fields[0]] = tuple(map(int, fields[1:]))
    return index


def fasta_extract_regions(fa_fname, intervals):
    """Extract an iterable of regions from an indexed FASTA file.

    Input: indexed FASTA file name; iterable of (seq_id, start, end) (1-based)
    Output: iterable of string sequences.
    """
    index = read_fasta_index(fa_fname)
    with open(fa_fname, 'rb') as fa_file:
        for chrom, rows in groupby(intervals, lambda cse: cse[0]):
            # Seek to chrom offset in FASTA
            try:
                seq_len, offset, chars_per_line, bytes_per_line = index[chrom]
            except KeyError:
                raise ValueError("Sequence ID '" + chrom + "' is not in FASTA "
                                 + "file " + fa_fname)
            echo("Extracting sequences from chromosome", chrom)
            eol_size = bytes_per_line - chars_per_line  # Handle \n\r, \n
            for _chrom, start, end in rows:
                if end > seq_len:
                    raise ValueError("Chromosome {} has length {} "
                                     .format(chrom, seq_len) +
                                     "which requested range {}:{}-{} "
                                     .format(_chrom, start, end) +
                                     "extends beyond")

                # Jump to the subsequence start position
                n_eols_to_skip, line_remainder = divmod(start, chars_per_line)
                skip_length = start + n_eols_to_skip * eol_size
                fa_file.seek(offset + skip_length)
                # Calculate how many bytes to read to capture the subsequence
                subseq_length = end - start
                line_to_go = chars_per_line - line_remainder
                n_eols_in_subseq = int(math.ceil((subseq_length - line_to_go)
                                                 / chars_per_line))
                # Read ahead by this many bytes
                subseq_bytes = fa_file.read(subseq_length
                                            + n_eols_in_subseq * eol_size)
                subseq = ''.join(subseq_bytes.split())  # Remove EOL characters
                # core.assert_equal("Number of characters read does not match "
                #                   "the number requested",
                #                   read=len(subseq),
                #                   requested=subseq_length)
                assert len(subseq) == subseq_length, (
                    "%s:%d-%d read bytes=%d, chars=%d; wanted chars=%d, eols=%d"
                    % (_chrom, start, end,
                       len(subseq_bytes), len(subseq),
                       subseq_length, n_eols_in_subseq))

                yield subseq


def _fasta_extract_regions_safe(fa_fname, intervals):
    """Simpler, slower version of fasta_extract_regions, for testing it."""
    from Bio import SeqIO
    idx = SeqIO.index(fa_fname, 'fasta')
    for chrom, rows in groupby(intervals, lambda cse: cse[0]):
        echo("Extracting sequences from chromosome", chrom)
        seq = str(idx[chrom].seq)
        for _chrom, start, end in rows:
            subseq = seq[start:end]
            if len(subseq) != end - start:
                echo("Short subsequence %s:%d-%d" % (_chrom, start, end),
                     "read", len(subseq), ", wanted", end - start)
            yield subseq
