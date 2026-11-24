"""Utilities."""

from __future__ import annotations

import contextlib
import functools
import gzip
import io
import os
from typing import IO, TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

#: Leading bytes of a gzip member (RFC 1952). BGZF, the block compression used
#: by tabix-indexed files, is a series of gzip members, so it begins with these
#: bytes too, and `gzip` reads all of its members in order.
GZIP_MAGIC = b"\x1f\x8b"


@contextlib.contextmanager
def as_text_handle(infile: str | os.PathLike[str] | IO[str]) -> Iterator[IO[str]]:
    """Open a path as a text handle, decompressing gzip transparently.

    Replaces `Bio.File.as_handle` for the readers that parse one line at a
    time. Compression is detected from the file's leading bytes rather than
    its name. The path is opened once and those bytes are peeked rather than
    consumed, so a named pipe or process substitution is still read whole.
    An uncompressed file is read exactly as ``open(path, "r")`` would read it.

    An already-open handle is yielded unchanged and left open; it belongs to
    the caller, who chose its mode.

    The readers built on `pandas.read_csv` do not need this, because pandas
    infers the compression from the file name, though only from suffixes it
    recognizes, such as ".gz" but not ".bgz".
    """
    if not isinstance(infile, str | os.PathLike):
        yield infile
        return
    with open(infile, "rb") as raw:
        if raw.peek(len(GZIP_MAGIC)).startswith(GZIP_MAGIC):
            with gzip.open(raw, "rt") as handle:
                yield handle
        else:
            with io.TextIOWrapper(raw) as handle:
                yield handle


def report_bad_line(line_parser: Callable) -> Callable:
    @functools.wraps(line_parser)
    def wrapper(line):
        try:
            return line_parser(line)
        except ValueError as exc:
            raise ValueError("Bad line: %r" % line) from exc

    return wrapper
