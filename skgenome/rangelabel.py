"""Handle text genomic ranges as named tuples.

A range specification should look like ``chromosome:start-end``, e.g.
``chr1:1234-5678``, with 1-indexed integer coordinates. We also allow
``chr1:1234-`` or ``chr1:-5678``, where missing start becomes 0 and missing end
becomes None.
"""

import collections
import re
from typing import Union
from collections.abc import Sequence

Region = collections.namedtuple("Region", "chromosome start end")
NamedRegion = collections.namedtuple("NamedRegion", "chromosome start end gene")

re_label = re.compile(r"(\w[\w.]*)?:(\d+)?-(\d+)?\s*(\S+)?")


def from_label(text: str, keep_gene: bool = True) -> Union[Region, NamedRegion]:
    """Parse a chromosomal range specification.

    Parameters
    ----------
    text : string
        Range specification, which should look like ``chr1:1234-5678`` or
        ``chr1:1234-`` or ``chr1:-5678``, where missing start becomes 0 and
        missing end becomes None.
    keep_gene : bool
        If True, include gene names as a 4th field where available; otherwise return a
        3-field Region of chromosomal coordinates without gene labels.
    """
    match = re_label.match(text)
    if not match:
        raise ValueError(
            f"Invalid range spec: {text} (should be like: chr1:2333000-2444000)"
        )

    chrom, start, end, gene = match.groups()
    start = int(start) - 1 if start else None
    end = int(end) if end else None
    if keep_gene:
        gene = gene or ""
        return NamedRegion(chrom, start, end, gene)
    return Region(chrom, start, end)


def to_label(row: Region) -> str:
    """Convert a Region tuple to a region label."""
    return f"{row.chromosome}:{row.start + 1}-{row.end}"


def unpack_range(a_range: Union[str, Sequence]) -> Region:
    """Extract chromosome, start, end from a string or tuple.

    Examples::

        "chr1" -> ("chr1", None, None)
        "chr1:100-123" -> ("chr1", 99, 123)
        ("chr1", 100, 123) -> ("chr1", 100, 123)
    """
    if not a_range:
        return Region(None, None, None)
    if isinstance(a_range, str):
        if ":" in a_range and "-" in a_range:
            return from_label(a_range, keep_gene=False)  # type: ignore
        return Region(a_range, None, None)
    if isinstance(a_range, (list, tuple)):
        if len(a_range) == 3:
            return Region(*a_range)
        if len(a_range) == 4:
            return Region(*a_range[:3])
    raise ValueError(f"Not a range: {a_range!r}")
