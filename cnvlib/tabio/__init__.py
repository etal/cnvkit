"""I/O for tabular formats of genomic data (regions or features).
"""
from __future__ import absolute_import, division, print_function

import logging
import sys

import pandas as pd
from Bio.File import as_handle

from .. import core, ngfrills
from ..gary import GenomicArray as GA
from ..cnary import CopyNumArray as CNA
from ..vary import VariantArray as VA
from . import bedio, picard, textcoord, vcfio


def read(infile, fmt="tab", into=None, sample_id=None, meta=None, **kwargs):
    """Read tabular data from a file or stream into a genome object.

    Supported formats:

    ======      ========    ====
    Format      Code        Into
    ------      --------    ----
    TSV         tab         Gary
    VCF         vcf         Vary
    BED         bed
    Int'l       interval    Gary
    ======      ========    ====

    :Parameters:
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
            Additional fields to add to metadata.

    """
    if meta is None:
        meta = {}
    meta.update(kwargs)
    if sample_id is None and "sample_id" not in meta:
        if isinstance(infile, basestring):
            sample_id = core.fbase(infile)
        else:
            sample_id = '<unknown>'
        meta["sample_id"] = sample_id
    if fmt not in READERS:
        raise ValueError("Unknown format: %s" % fmt)
    try:
        reader, suggest_into = READERS[fmt]
        dframe = reader(infile)
    except ValueError:
        # File is blank/empty, most likely
        logging.info("Blank file?: %s", infile)
        suggest_into = GA
        dframe = []
    # TODO/ENH CategoricalIndex ---
    # if dframe:
    # dframe['chromosome'] = pd.Categorical(dframe['chromosome'],
    #                                      dframe.chromosome.drop_duplicates(),
    #                                      ordered=True)
    # Create a multi-index of genomic coordinates (like GRanges)
    # dframe.set_index(['chromosome', 'start'], inplace=True)
    return (into or suggest_into)(dframe, meta)


def read_auto(infile):
    """Auto-detect file format, and return an appropriate parser function."""
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


def read_cna(infile, sample_id=None, meta=None, **kwargs):
    """Read a tabular file to create a CopyNumArray object."""
    return read(infile, into=CNA, sample_id=sample_id, meta=meta, **kwargs)


def read_tab(infile):
    """Read tab-separated data with column names in the first row.

    The format is BED-like, but with a header row included and with
    arbitrary extra columns.
    """
    dframe = pd.read_table(infile, dtype={'chromosome': 'str'})
    if "log2" in dframe.columns:
        # Every bin needs a log2 value; the others can be NaN
        d2 = dframe.dropna(subset=["log2"])
        if len(d2) < len(dframe):
            logging.warn("Dropped %d rows with missing log2 values",
                        len(dframe) - len(d2))
            dframe = d2
    return dframe


READERS = {
    "auto": (read_auto, GA),
    "tab": (read_tab, GA),
    "bed": (bedio.read_bed, GA),
    "bed3": (bedio.read_bed3, GA),
    "bed4": (bedio.read_bed4, GA),
    "bed6": (bedio.read_bed6, GA),
    "interval": (picard.read_interval, GA),
    "text": (textcoord.read_text, GA),
    "vcf": (vcfio.read_vcf, VA),
}


# _____________________________________________________________________

def write(garr, outfile=None, fmt="tab", verbose=True):
    """Write a genome object to a file or stream."""
    dframe = WRITERS[fmt](garr.data)
    with ngfrills.safe_write(outfile or sys.stdout) as handle:
        dframe.to_csv(handle, header=(fmt == "tab"), index=False, sep='\t',
                      float_format='%.6g')
    if verbose:
        # Log the output path, if possible
        if isinstance(outfile, basestring):
            outfname = outfile
        elif hasattr(outfile, 'name') and outfile not in (sys.stdout,
                                                            sys.stderr):
            outfname = outfile.gene
        else:
            # Probably stdout or stderr used in a pipeline -- don't pollute
            return
        logging.info("Wrote %s with %d regions", outfname, len(dframe))


def write_tab(dframe):
    """Write tab-separated data with column names in the first row."""
    return dframe


WRITERS = {
    "tab": write_tab,
    "bed": bedio.write_bed,
    "interval": picard.write_interval,
    "text": textcoord.write_text,
    "vcf": vcfio.write_vcf,
}
