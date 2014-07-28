"""NGS utilities."""
from __future__ import division, print_function
import collections
import contextlib
import functools
import math
import os
import re
import shlex
import subprocess
import sys
import tempfile
from itertools import groupby, islice

from Bio.File import as_handle
from Bio._py3k import basestring, map

import pysam

def echo(*words):
    print(*words, file=sys.stderr)


# __________________________________________________________
# BED/Interval I/O

def sniff_num_columns(bed_fname):
    """Sniff the number of columns in a BED/interval file.

    Guidance:
        3 cols => coordinates only;
        5 cols => intervals file (coordinates, strand, name);
        otherwise => Full or extended BED format
    """
    for firstrow in parse_regions(bed_fname):
        return len(firstrow)


def sniff_region_format(fname):
    """Guess whether the file format is BED, Picard interval list, or text.

    Returns a tuple of the format name (str) and a handle to the opened file.
    """
    with open(fname, 'rU') as handle:
        line = handle.readline()
    if '\t' not in line and ':' in line and '-' in line:
        return 'text'
    if line.startswith('@') or re.match('\w+\t\d+\t\d+\t(\+|-|\.)\t\w+', line):
        echo("Sniffed interval")
        return 'interval'
    if line.startswith('track') or line.count('\t') > 1:
        echo("Sniffed BED")
        return 'bed'
    raise ValueError("WTF format is this?: %s" % line)


def parse_regions(fname, coord_only=False, keep_strand=False):
    """Parse regions in any of the expected file formats.

    Iterates over tuples of the tabular contents. Header lines are skipped.

    Start and end coordinates are base-0, half-open.

    If coord_only, yield triplets of (chrom, start, end). Otherwise, yield
    quads of (chrom, start, end, name).
    """
    fmt = sniff_region_format(fname)
    parser = {'text': parse_text_coords,
              'interval': parse_interval_list,
              'bed': parse_bed,
             }[fmt]
    return parser(fname, coord_only, keep_strand)


def report_bad_line(line_parser):
    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError:
            raise ValueError("Bad line: %r" % line)
    return wrapper


def parse_text_coords(fname, coord_only, keep_strand):
    """Parse text coordinates: chrom:start-end

    Text coordinates are assumed to be counting from 1.
    """
    if coord_only:
        @report_bad_line
        def _parse_line(line):
            chrom, _rest = line.rstrip().split(':', 1)
            start, end = _rest.split('-')
            if ':' in end:
                end = end.split(':', 1)[0]
            return chrom, int(start) - 1, int(end)
    else:
        @report_bad_line
        def _parse_line(line):
            fields = line.split(':')
            if len(fields) == 3:
                chrom, start_end, name = fields
            elif len(fields) == 2:
                chrom, start_end = fields
                name = ''
            else:
                raise ValueError
            start, end = start_end.split('-')
            return chrom, int(start) - 1, int(end), name.rstrip()

    with as_handle(fname, 'rU') as handle:
        for line in handle:
            yield _parse_line(line)


def parse_interval_list(fname, coord_only, keep_strand):
    """Parse a Picard-compatible interval list.

    Expected tabular columns:
        chromosome, start position, end position, strand, region name

    Counting is from 1.
    """
    if coord_only:
        if keep_strand:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end, strand = line.split('\t')[:4]
                return chrom, int(start) - 1, int(end), strand.rstrip()
        else:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end = line.split('\t')[:3]
                return chrom, int(start) - 1, int(end)
    elif keep_strand:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t')
            chrom, start, end, strand = fields[:4]
            if len(fields) > 4:
                name = fields[-1].rstrip()
            else:
                name = ''
            return chrom, int(start) - 1, int(end), name, strand
    else:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t')
            chrom, start, end = fields[:3]
            if len(fields) > 3:
                name = fields[-1].rstrip()
            else:
                name = ''
            return chrom, int(start) - 1, int(end), name

    with as_handle(fname, 'rU') as handle:
        for line in handle:
            if line.startswith('@'):
                # Skip the SAM header
                continue
            yield _parse_line(line)


def parse_bed(fname, coord_only, keep_strand):
    """Parse a BED file.

    A BED file has these columns:
        chromosome, start position, end position, [name, strand, other stuff...]

    Counting is from 0.

    Sets of regions are separated by "track" lines. This function stops
    iteration after encountering a track line other than the first one in the
    file.
    """
    if coord_only:
        if keep_strand:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end, _name, _score, strand = line.split('\t', 6)[:6]
                return chrom, int(start), int(end), strand.rstrip()
        else:
            @report_bad_line
            def _parse_line(line):
                chrom, start, end = line.split('\t', 3)[:3]
                return chrom, int(start), int(end)
    elif keep_strand:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t', 6)
            chrom, start, end = fields[:3]
            name = (fields[3].rstrip()
                    if len(fields) >= 4 else '')
            strand = (fields[5].rstrip()
                      if len(fields) >= 6 else '.')
            return chrom, int(start), int(end), name, strand
    else:
        @report_bad_line
        def _parse_line(line):
            fields = line.split('\t', 4)
            chrom, start, end = fields[:3]
            name = (fields[3].rstrip()
                    if len(fields) >= 4 else '')
            return chrom, int(start), int(end), name

    with as_handle(fname, 'rU') as handle:
        firstline = next(handle)
        if firstline.startswith("track"):
            pass
        else:
            yield _parse_line(firstline)

        for line in handle:
            if line.startswith('track'):
                raise StopIteration
            yield _parse_line(line)


def parse_bed_track(line):
    """Parse the "name" field of a BED track definition line.

    Example:
    track name=146793_BastianLabv2_P2_target_region description="146793_BastianLabv2_P2_target_region"
    """
    fields = shlex.split(line)  # raises ValueError if line is corrupted
    assert fields[0] == 'track'
    for field in fields[1:]:
        if '=' in field:
            key, val = field.split('=', 1)
            if key == 'name':
                return val
    else:
        raise ValueError("No name defined for this track")


def group_bed_tracks(bedfile):
    """Group the parsed rows in a BED file by track.

    Yields (track_name, iterable_of_lines), much like itertools.groupby.
    """
    # ENH - make this memory-efficient w/ generators or something
    with as_handle(bedfile, 'r') as handle:
        curr_track = 'DEFAULT'
        curr_lines = []
        for line in handle:
            if line.startswith('track'):
                if curr_lines:
                    yield curr_track, curr_lines
                    curr_lines = []
                curr_track = parse_bed_track(line)
            else:
                curr_lines.append(line)
        yield curr_track, curr_lines


# __________________________________________________________
# Indexed FASTA I/O

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
    # Build a dict of keys -> offsets
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
                _seq_len, offset, chars_per_line, bytes_per_line = index[chrom]
            except KeyError:
                raise ValueError("Sequence ID '" + chrom + "' is not in FASTA "
                                 + "file " + fa_fname)
            eol_size = bytes_per_line - chars_per_line  # Handle \n\r, \n
            for _chrom, start, end in rows:
                start -= 1
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
                    "Read bytes=%d, chars=%d; wanted chars=%d, eols=%d"
                    % (len(subseq_bytes), len(subseq),
                       subseq_length, n_eols_in_subseq))

                yield subseq


def _fasta_extract_regions_safe(fa_fname, intervals):
    """Simpler, slower version of fasta_extract_regions, for testing it."""
    from Bio import SeqIO
    idx = SeqIO.index(fa_fname, 'fasta')
    for chrom, rows in groupby(intervals, lambda cse: cse[0]):
        seq = str(idx[chrom].seq)
        for _chrom, start, end in rows:
            start -= 1
            yield seq[start:end]


# __________________________________________________________
# VCF I/O

def load_vcf(fname, min_depth=1, skip_hom=True):
    """Parse SNV coordinates from a VCF file; group by chromosome.

    Returns a dict of: {chrom: position, zygosity, alt allele frequency}
    """
    chrom_snvs = collections.defaultdict(list)
    for record, zygosity in filter_vcf_lines(fname, min_depth, skip_hom):
        # Read depth
        samp = record.samples[0]
        depth = float(samp.data.DP)
        alt_count = float(samp.data.AD[1])
        alt_freq = alt_count / depth
        chrom_snvs[record.CHROM].append((record.POS, zygosity, alt_freq))
    return chrom_snvs


def filter_vcf_lines(vcf_fname, min_depth, skip_hom):
    import vcf
    with open(vcf_fname) as vcffile:
        vcf_reader = vcf.Reader(vcffile)
        for record in vcf_reader:
            # Skip unassigned contigs
            if len(record.CHROM) > len("chr22"):
                continue
            # Skip homozygous variants (optionally)
            # XXX if no 'AF', check 'FA' (for MuTect)
            zygosity = record.INFO['AF'][0]  # 1.0 or 0.5
            if skip_hom and zygosity != 0.5:
                continue
            samp = record.samples[0]
            depth = samp.data.DP
            if depth < min_depth:
                continue
            yield record, zygosity



# __________________________________________________________
# Shell

def call_quiet(*args):
    """Safely run a command and get stdout; print stderr if there's an error.

    Like subprocess.check_output, but silent in the normal case where the
    command logs unimportant stuff to stderr. If there is an error, then the
    full error message(s) is shown in the exception message.
    """
    # args = map(str, args)
    if not len(args):
        raise ValueError("Must supply at least one argument (the command name)")
    try:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except OSError, exc:
        raise RuntimeError("Could not find the executable %r" % args[0]
                           + " -- is it installed correctly?"
                           + "\n(Original error: %s)" % exc)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("Subprocess command failed:\n$ %s\n\n%s"
                           % (' '.join(args), err))
    return out


def ensure_bam_index(bam_fname):
    """Ensure a BAM file is indexed, to enable fast traversal & lookup."""
    bai_fname = bam_fname + '.bai'
    if not is_newer_than(bai_fname, bam_fname):
        echo("Indexing BAM file", bam_fname)
        pysam.index(bam_fname)
    assert os.path.isfile(bai_fname), "Failed to generate index " + bai_fname
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
        getter = lambda read: read.qname
    else:
        # Compare read locations
        getter = lambda read: read.pos

    # TODO - repeat at 50%, ~99% through the BAM
    bam = pysam.Samfile(bam_fname, 'rb')
    last_read = None
    for read in islice(bam, span):
        read_value = getter(read)
        if last_read is not None:
            if read_value < last_read:
                return False
        last_read = read_value

    return True


def ensure_fasta_index(fasta_fname):
    """Ensure a FASTA file is indexed for samtools, to enable fast lookup."""
    fai_fname = fasta_fname + '.fai'
    if not is_newer_than(fai_fname, fasta_fname):
        echo("Indexing FASTA file", fasta_fname)
        pysam.faidx(fasta_fname)
    assert os.path.isfile(fai_fname), "Failed to generate index " + fai_fname
    return fai_fname


def is_newer_than(target_fname, orig_fname):
    if not os.path.isfile(target_fname):
        return False
    return (os.stat(target_fname).st_mtime >= os.stat(orig_fname).st_mtime)


def ensure_path(fname):
    """Create dirs and move an existing file to avoid overwriting, if necessary.

    If a file already exists at the given path, it is renamed with an integer
    suffix to clear the way.
    """
    if '/' in os.path.normpath(fname):
        # Ensure the output directory exists
        dname = os.path.dirname(os.path.abspath(fname))
        if dname and not os.path.isdir(dname):
            try:
                os.makedirs(dname)
            except OSError, exc:
                raise OSError("Output path " + fname +
                                " contains a directory " + dname +
                                " that cannot be created: " + str(exc))
    if os.path.isfile(fname):
        # Add an integer suffix to the existing file name
        cnt = 1
        bak_fname = "%s.%d" % (fname, cnt)
        while os.path.isfile(bak_fname):
            cnt += 1
            bak_fname = "%s.%d" % (fname, cnt)
        os.rename(fname, bak_fname)
        echo("Moved existing file", fname, "->", bak_fname)
    return True


@contextlib.contextmanager
def safe_write(outfile, verbose=True):
    """Write to a filename or file-like object with error handling.

    If given a file name, open it. If the path includes directories that don't
    exist yet, create them.  If given a file-like object, just pass it through.
    """
    if isinstance(outfile, basestring):
        dirname = os.path.dirname(outfile)
        if dirname and not os.path.isdir(dirname):
            os.mkdir(dirname)
            echo("Created directory", dirname)
        with open(outfile, 'w') as handle:
            yield handle
    else:
        yield outfile

    # Log the output path, if possible
    if verbose:
        if isinstance(outfile, basestring):
            outfname = outfile
        elif hasattr(outfile, 'name') and outfile not in (sys.stdout,
                                                          sys.stderr):
            outfname = outfile.name
        else:
            # Probably stdout or stderr -- don't ruin the pipeline
            return
        echo("Wrote", outfname)


@contextlib.contextmanager
def temp_write_text(text):
    """Save text to a temporary file.

    NB: This won't work on Windows b/c the file stays open.
    """
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write(text)
        tmp.flush()
        yield tmp.name
