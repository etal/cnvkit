"""I/O for tabular formats of genomic data (regions or features).
"""
from __future__ import absolute_import, division, print_function
from past.builtins import basestring

import logging
import sys

import pandas as pd
from Bio.File import as_handle

from .. import core, ngfrills
from ..gary import GenomicArray as GA
from ..cnary import CopyNumArray as CNA
from ..vary import VariantArray as VA
from . import bedio, picard, seg, tab, textcoord, vcfio


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
            meta["sample_id"] = core.fbase(infile)
        elif hasattr(infile, "name"):
            meta["sample_id"] = core.fbase(infile.name)
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

    fmt = ngfrills.sniff_region_format(infile)
    if hasattr(infile, "seek"):
        infile.seek(0)
    if fmt is None:
        fmt = "tab"
    if fmt == "bed":
        logging.info("Detected file format: BED")
    elif fmt == "interval":
        logging.info("Detected file format: interval list")
    return read(infile, fmt)


def read_cna(infile, sample_id=None, meta=None):
    """Read a tabular file to create a CopyNumArray object."""
    return read(infile, into=CNA, sample_id=sample_id, meta=meta)


READERS = {
    # Format name, formatter, default target class
    "auto": (read_auto, GA),
    "bed": (bedio.read_bed, GA),
    "bed3": (bedio.read_bed3, GA),
    "bed4": (bedio.read_bed4, GA),
    "bed6": (bedio.read_bed6, GA),
    "interval": (picard.read_interval, GA),
    "picardhs": (picard.read_picard_hs, CNA),
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
    with core.safe_write(outfile or sys.stdout, verbose=False) as handle:
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
    "interval": (picard.write_interval, False),
    "picardhs": (picard.write_picard_hs, True),
    "seg": (seg.write_seg, True),
    "tab": (tab.write_tab, True),
    "text": (textcoord.write_text, False),
    "vcf": (vcfio.write_vcf, True),
}
