"""Utilities."""

from __future__ import annotations
import functools
from typing import TYPE_CHECKING
from collections.abc import Callable


def report_bad_line(line_parser: Callable) -> Callable:
    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError as exc:
            raise ValueError("Bad line: %r" % line) from exc

    return wrapper
