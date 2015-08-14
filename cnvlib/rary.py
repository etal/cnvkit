"""An array of genomic regions or features."""
from __future__ import absolute_import, division, print_function

import pandas as pd
from Bio.File import as_handle

from . import core, gary
from .ngfrills import echo
from .ngfrills.regions import sniff_region_format, report_bad_line


class RegionArray(gary.GenomicArray):
    """An array of genomic intervals."""
    _required_columns = ("chromosome", "start", "end",
                         # "name", "strand",
                        )

    def __init__(self, data_table, meta_dict=None):
        gary.GenomicArray.__init__(self, data_table, meta_dict)

    @classmethod
    def read(cls, fname, sample_id=None, fmt=None):
        """Read regions in any of the expected file formats.

        Iterates over tuples of the tabular contents. Header lines are skipped.

        Start and end coordinates are base-0, half-open.
        """
        if sample_id is None:
            if isinstance(fname, basestring):
                sample_id = core.fbase(fname)
            elif fmt is None:
                raise ValueError("To read regions from a stream, the file "
                                 "format must be specified with the `fmt` "
                                 "argument.")
            else:
                sample_id = '<unknown>'

        if not fmt:
            fmt = sniff_region_format(fname)
            if fmt is None:
                return cls([])
            if fmt == 'bed':
                echo("Detected file format: BED")
            elif fmt == 'interval':
                echo("Detected file format: interval list")
        parser = {'text': _parse_text_coords,
                  'interval': _parse_interval_list,
                  'bed': _parse_bed,
                 }[fmt]
        table = parser(fname)
        return cls(table, {"sample_id": sample_id})


def _parse_text_coords(infile):
    """Parse text coordinates: chrom:start-end

    Or sometimes: chrom:start-end:name

    Text coordinates are assumed to be counting from 1.
    """
    @report_bad_line
    def _parse_line(line):
        fields = line.split(':')
        if len(fields) == 3:
            chrom, start_end, name = fields
        elif len(fields) == 2:
            chrom, start_end = fields
            name = '-'
        else:
            raise ValueError("Bad line: %r" % line)
        start, end = start_end.split('-')
        return chrom, int(start) - 1, int(end), name.rstrip()

    with as_handle(infile, 'rU') as handle:
        rows = [_parse_line(line) for line in handle]
    return pd.DataFrame.from_records(rows, columns=["chromosome", "start",
                                                    "end", "name"])


def _parse_interval_list(infile):
    """Parse a Picard-compatible interval list.

    Expected tabular columns:
        chromosome, start position, end position, strand, region name

    Counting is from 1.
    """
    table = pd.read_table(infile,
                          comment='@', # Skip the SAM header
                          names=["chromosome", "start", "end", "name",
                                 "strand"])
    table["name"].fillna('-', inplace=True)
    table["start"] -= 1
    return table


def _parse_bed(infile):
    """Parse a BED file.

    A BED file has these columns:
        chromosome, start position, end position, [name, strand, other stuff...]

    Counting is from 0.

    Sets of regions are separated by "track" lines. This function stops reading
    after encountering a track line other than the first one in the file.
    """
    # ENH: just pd.read_table, skip 'track'
    @report_bad_line
    def _parse_line(line):
        fields = line.split('\t', 6)
        chrom, start, end = fields[:3]
        name = (fields[3].rstrip()
                if len(fields) >= 4 else '-')
        strand = (fields[5].rstrip()
                if len(fields) >= 6 else '.')
        return chrom, int(start), int(end), name, strand

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
                                                    "end", "name", "strand"])

