"""I/O for tabular formats of genomic data.
"""
from __future__ import print_function, absolute_import, division

import logging
import sys

import pandas as pd
from Bio.File import as_handle

from cnvlib import core, ngfrills
from cnvlib.gary import GenomicArray as GA
from cnvlib.vary import VariantArray as VA


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


def read_sniff(infile):
    """Auto-detect file format, and return an appropriate parser function."""
    if not isinstance(infile, basestring) and not hasattr(infile, "seek"):
        raise ValueError("Can only auto-detect format from filename or "
                         "seekable (local, on-disk) files, which %s is not"
                         % infile)

    fmt = ngfrills.sniff_region_format(infile)
    if fmt is None:
        return None #cls([])
    if fmt == 'bed':
        logging.info("Detected file format: BED")
    elif fmt == 'interval':
        logging.info("Detected file format: interval list")
    parser = {'text': read_text,
              'interval': read_interval,
              'bed': read_bed,
             }[fmt]
    return parser


def read_tab(infile):
    """Read tab-separated data with column names in the first row."""
    try:
        dframe = pd.read_table(infile, #na_filter=False,
                            dtype={'chromosome': 'string'})
        if "log2" in dframe.columns:
            # Every bin needs a log2 value; the others can be NaN
            t2 = dframe.dropna(subset=["log2"])
            if len(t2) < len(dframe):
                logging.warn("Dropped %d rows with missing log2 values",
                            len(dframe) - len(t2))
                dframe = t2
    except ValueError:
        # File is blank/empty, most likely
        logging.warn("Blank file?: %s", infile)
        dframe = []
    return dframe


def read_bed(infile):
    """UCSC Browser Extensible Data (BED) format.

    A BED file has these columns:
        chromosome, start position, end position, [gene, strand, other stuff...]

    Coordinate indexing is from 0.

    Sets of regions are separated by "track" lines. This function stops reading
    after encountering a track line other than the first one in the file.
    """
    # ENH: just pd.read_table, skip 'track'
    @ngfrills.report_bad_line
    def _parse_line(line):
        fields = line.split('\t', 6)
        chrom, start, end = fields[:3]
        gene = (fields[3].rstrip()
                if len(fields) >= 4 else '-')
        strand = (fields[5].rstrip()
                if len(fields) >= 6 else '.')
        return chrom, int(start), int(end), gene, strand

    def track2track(handle):
        firstline = next(handle)
        if firstline.startswith("track"):
            pass
        else:
            yield firstline
        for line in handle:
            if line.startswith('track'):
                raise StopIteration
            yield line

    with as_handle(infile, 'rU') as handle:
        rows = map(_parse_line, track2track(handle))
    return pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                    "end", "gene", "strand"])


# TODO - would these be useful?
def read_bed3(infile):
    """3-column BED format: chromosome, start, end."""
    return NotImplemented

def read_bed4(infile):
    """4-column BED format: chromosome, start, end, name."""
    return NotImplemented

def read_bed6(infile):
    """6-column BED format: chromosome, start, end, name, score, strand."""
    return NotImplemented


def read_interval(infile):
    """GATK/Picard-compatible interval list format.

    Expected tabular columns:
        chromosome, start position, end position, strand, gene

    Coordinate indexing is from 1.
    """
    dframe = pd.read_table(infile,
                           comment='@', # Skip the SAM header
                           names=["chromosome", "start", "end", "strand", "gene",
                                 ])
    dframe["gene"].fillna('-', inplace=True)
    dframe["start"] -= 1
    return dframe


def read_text(infile):
    """Text coordinate format: "chr:start-end", one per line.

    Or sometimes: "chrom:start-end gene" or "chrom:start-end REF>ALT"

    Coordinate indexing is assumed to be from 1.
    """
    @ngfrills.report_bad_line
    def _parse_line(line):
        fields = line.split(':')
        if len(fields) == 3:
            chrom, start_end, gene = fields
        elif len(fields) == 2:
            chrom, start_end = fields
            gene = '-'
        else:
            raise ValueError("Bad line: %r" % line)
        start, end = start_end.split('-')
        return chrom, int(start) - 1, int(end), gene.rstrip()

    with as_handle(infile, 'rU') as handle:
        rows = [_parse_line(line) for line in handle]
    return pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                    "end", "gene"])


def read_vcf(infile):
    """Variant Call Format (VCF) for SNV loci."""
    return NotImplemented


READERS = {
    "tab": (read_tab, GA),
    "bed": (read_bed, GA),
    "bed3": (read_bed3, GA),
    "bed4": (read_bed4, GA),
    "bed6": (read_bed6, GA),
    "interval": (read_interval, GA),
    # "interval_list": (read_interval, GA),
    "sniff": (read_sniff, GA),
    "text": (read_text, GA),
    # "text_coords": (read_text, GA),
    "vcf": (read_vcf, VA),
}



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


def write_bed(dframe):
    if len(dframe.columns) == 3:
        return write_bed3(dframe)
    elif len(dframe.columns) == 3:
        return write_bed4(dframe)
    else:
        # Default: BED-like, keep all trailing columns
        return dframe


def write_bed3(dframe):
    return dframe.loc[:, ["chromosome", "start", "end"]]


def write_bed4(dframe):
    dframe = dframe.copy()
    if "gene" not in dframe:
        dframe["gene"] = '-'
    return dframe.loc[:, ["chromosome", "start", "end", "gene"]]


def write_interval(dframe):
    dframe = dframe.copy()
    dframe["start"] += 1
    if "gene" not in dframe:
        dframe["gene"] = '-'
    if "strand" not in dframe:
        dframe["strand"] = "+"
    return dframe.loc[:, ["chromosome", "start", "end", "strand", "gene"]]


def write_text(dframe):
    dframe = dframe.copy()
    dframe['start'] += 1
    return dframe.apply(GA.row2label, axis=1)


def write_vcf(infile):
    """Variant Call Format (VCF) for SV loci."""
    return NotImplemented


WRITERS = {
    "tab": write_tab,
    "bed": write_bed,
    "interval": write_interval,
    "text": write_text,
    "vcf": write_vcf,
}

