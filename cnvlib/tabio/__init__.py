"""I/O for tabular formats of genomic data (regions or features).
"""
from __future__ import absolute_import, division, print_function
from past.builtins import basestring

import collections
import logging
import re
import sys

import pandas as pd
from Bio.File import as_handle

from skgenome import GenomicArray as GA
from . import bedio, genepred, gff, picard, seg, tab, textcoord, vcfio

from ..core import fbase, safe_write
from ..cnary import CopyNumArray as CNA
from ..vary import VariantArray as VA


def read(infile, fmt="tab", into=None, sample_id=None, meta=None, **kwargs):
    """Read tabular data from a file or stream into a genome object.

    Supported formats: see `READERS`

    If a format supports multiple samples, return the sample specified by
    `sample_id`, or if unspecified, return the first sample and warn if there
    were other samples present in the file.

    Parameters
    ----------
    infile : handle or string
        Filename or opened file-like object to read.
    fmt : string
        File format.
    into : class
        GenomicArray class or subclass to instantiate, overriding the
        default for the target file format.
    sample_id : string
        Sample identifier.
    meta : dict
        Metadata, as arbitrary key-value pairs.
    **kwargs :
        Additional keyword arguments to the format-specific reader function.

    Returns
    -------
    GenomicArray or subclass
        The data from the given file instantiated as `into`, if specified, or
        the default base class for the given file format (usually GenomicArray).
    """
    if fmt == 'auto':
        return read_auto(infile)
    elif fmt in READERS:
        reader, suggest_into = READERS[fmt]
    else:
        raise ValueError("Unknown format: %s" % fmt)

    if meta is None:
        meta = {}
    if "sample_id" not in meta:
        if sample_id:
            meta["sample_id"] = sample_id
        elif isinstance(infile, basestring):
            meta["sample_id"] = fbase(infile)
        elif hasattr(infile, "name"):
            meta["sample_id"] = fbase(infile.name)
        else:
            # meta["sample_id"] = "<unknown>"
            pass
    if "filename" not in meta:
        if isinstance(infile, basestring):
            meta["filename"] = infile
        elif hasattr(infile, "name"):
            meta["filename"] = infile.name
    if fmt in ("seg", "vcf") and sample_id is not None:
        # Multi-sample formats: choose one sample
        kwargs["sample_id"] = sample_id
    try:
        dframe = reader(infile, **kwargs)
    except pd.io.common.EmptyDataError:
        # File is blank/empty, most likely
        logging.info("Blank %s file?: %s", fmt, infile)
        dframe = []
    result = (into or suggest_into)(dframe, meta)
    result.sort_columns()
    result.sort()
    return result
    # ENH CategoricalIndex ---
    # if dframe:
    # dframe['chromosome'] = pd.Categorical(dframe['chromosome'],
    #                                      dframe.chromosome.drop_duplicates(),
    #                                      ordered=True)
    # Create a multi-index of genomic coordinates (like GRanges)
    # dframe.set_index(['chromosome', 'start'], inplace=True)


def read_auto(infile):
    """Auto-detect a file's format and use an appropriate parser to read it."""
    if not isinstance(infile, basestring) and not hasattr(infile, "seek"):
        raise ValueError("Can only auto-detect format from filename or "
                         "seekable (local, on-disk) files, which %s is not"
                         % infile)

    fmt = sniff_region_format(infile)
    if hasattr(infile, "seek"):
        infile.seek(0)
    if fmt:
        logging.info("Detected file format: " + fmt)
    else:
        # File is blank -- simple BED will handle this OK
        fmt = "bed3"
    return read(infile, fmt or 'tab')


READERS = {
    # Format name, formatter, default target class
    "auto": (read_auto, GA),
    "bed": (bedio.read_bed, GA),
    "bed3": (bedio.read_bed3, GA),
    "bed4": (bedio.read_bed4, GA),
    "bed6": (bedio.read_bed6, GA),
    "gff": (gff.read_gff, GA),
    "interval": (picard.read_interval, GA),
    "refflat": (genepred.read_refflat, GA),
    "picardhs": (picard.read_picard_hs, GA),
    "seg": (seg.read_seg, CNA),
    "tab": (tab.read_tab, GA),
    "text": (textcoord.read_text, GA),
    "vcf": (vcfio.read_vcf, VA),
}


# _____________________________________________________________________

def write(garr, outfile=None, fmt="tab", verbose=True, **kwargs):
    """Write a genome object to a file or stream."""
    formatter, show_header = WRITERS[fmt]
    if fmt in ("seg", "vcf"):
        kwargs["sample_id"] = garr.sample_id
    dframe = formatter(garr.data, **kwargs)
    with safe_write(outfile or sys.stdout, verbose=False) as handle:
        dframe.to_csv(handle, header=show_header, index=False, sep='\t',
                      float_format='%.6g')
    if verbose:
        # Log the output path, if possible
        if isinstance(outfile, basestring):
            outfname = outfile
        elif hasattr(outfile, 'name') and outfile not in (sys.stdout,
                                                          sys.stderr):
            outfname = outfile.name
        else:
            # Probably stdout or stderr used in a pipeline -- don't pollute
            return
        logging.info("Wrote %s with %d regions", outfname, len(dframe))


WRITERS = {
    # Format name, formatter, show header
    "bed": (bedio.write_bed, False),
    "bed3": (bedio.write_bed3, False),
    "bed4": (bedio.write_bed4, False),
    # "gff": (gff.write_gff, False),
    "interval": (picard.write_interval, False),
    "picardhs": (picard.write_picard_hs, True),
    "seg": (seg.write_seg, True),
    "tab": (tab.write_tab, True),
    "text": (textcoord.write_text, False),
    "vcf": (vcfio.write_vcf, True),
}


# _____________________________________________________________________

def sniff_region_format(fname):
    """Guess the format of the given file by reading the first line.

    Returns
    -------
    str or None
        The detected format name, or None if the file is empty.
    """
    # has_track = False
    with as_handle(fname, 'rU') as handle:
        for line in handle:
            if not line.strip():
                # Skip blank lines
                continue
            if line.startswith('track'):
                # NB: Could be UCSC BED or Ensembl GFF
                # has_track = True
                continue
            # Formats that (may) declare themselves in an initial '#' comment
            if (line.startswith('##gff-version') or
                format_patterns['gff'].match(line)):
                return 'gff'
            if line.startswith(('##fileformat=VCF', '#CHROM\tPOS\tID')):
                return 'vcf'
            if line.startswith('#'):
                continue
            # Formats that need to be guessed solely by regex
            if format_patterns['text'].match(line):
                return 'text'
            if format_patterns['tab'].match(line):
                return 'tab'
            if line.startswith('@') or format_patterns['interval'].match(line):
                return 'interval'
            if format_patterns['refflat'].match(line):
                return 'refflat'
            if format_patterns['bed'].match(line):
                return 'bed'

            raise ValueError("File %r does not appear to be a recognized "
                             "format! (Any of: %s)\n"
                             "First non-blank line:\n%s"
                             % (fname, ', '.join(format_patterns.keys()), line))


format_patterns = collections.OrderedDict([
    #  ('genepred', re.compile()),
    #  ('genepredext', re.compile()),
    ('text', re.compile(r'\w+:\d*-\d*.*')),
    ('tab', re.compile('\t'.join(('chromosome', 'start', 'end')))),
    ('interval', re.compile('\t'.join((
        r'\w+', r'\d+', r'\d+', r'[.+-]', r'\S+$')))),
    ('refflat', re.compile('\t'.join((
        r'\S+', r'\S+', r'\w+', r'[+-]',
        r'\d+', r'\d+', r'\d+', r'\d+', r'\d+',
        r'(\d+,)+', r'(\d+,)+$')))),
    ('gff', re.compile('\t'.join((
        r'\w+', r'\S+', r'\w+', r'\d+', r'\d+',
        r'\S+', r'[.?+-]', r'[012.]', r'.*')))),
    ('bed', re.compile('\t'.join((r'\w+', r'\d+', r'\d+')))),
])
