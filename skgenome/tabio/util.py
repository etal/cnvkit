"""Utilities."""
import functools


def report_bad_line(line_parser):
    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError as exc:
            raise ValueError("Bad line: %r" % line) from exc

    return wrapper
