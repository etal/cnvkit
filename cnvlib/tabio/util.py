"""Utilities."""
from __future__ import absolute_import, division, print_function

import functools


def report_bad_line(line_parser):
    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError:
            raise ValueError("Bad line: %r" % line)
    return wrapper
