"""Utilities."""
from __future__ import absolute_import, division, print_function

import functools
import re

from Bio.File import as_handle


def report_bad_line(line_parser):
    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError:
            raise ValueError("Bad line: %r" % line)
    return wrapper


def sniff_region_format(fname):
    """Guess whether the file format is BED, Picard interval list, or text.

    Returns a tuple of the format name (str) or None if the file is empty.
    """
    with as_handle(fname, 'rU') as handle:
        for line in handle:
            if not line.strip():
                continue
            if '\t' not in line and ':' in line and '-' in line:
                return 'text'
            if line.startswith('@') or re.match('\w+\t\d+\t\d+\t(\+|-|\.)\t\S+',
                                                line):
                return 'interval'
            if line.startswith('track') or line.count('\t') > 1:
                return 'bed'
            raise ValueError("File " + repr(fname) + " does not appear to "
                             + "be BED, interval list, or 'chr:start-end' "
                             + "text!\nFirst non-blank line: " + repr(line))
