"""How far into the reference a VCF record reaches.

A record's span is not simply its POS and REF. The specification defines it as
the maximum of the reference allele's length and the lengths implied by
INFO/END and by INFO/SVLEN on symbolic alleles. htslib computes that maximum
from 1.23 onward; through 1.21 -- the version pysam 0.23.3 bundles, which is
this project's pinned floor -- it honoured a declared END even where the
reference allele reached further, and never consulted SVLEN at all.

CNVkit therefore takes the maximum itself rather than inheriting whichever
answer the installed pysam happens to give, so that a run's coordinates do not
depend on which wheel is on the system. Only the SVLEN half needs sharing
between the readers: `vcfio` can take htslib's own `record.stop` as a floor
and raise it, while the text readers in `vcfsimple` build the span from the
raw fields, but the two must agree on which alleles let SVLEN reach the
reference at all.

The gVCF FORMAT/LEN term is deliberately absent. It sizes `<*>` and
`<NON_REF>` reference blocks, which `vcfio` drops as non-alleles and the text
readers cannot see, since they read no sample columns; where a newer htslib
has applied it, taking `record.stop` as a floor preserves it.

One case reaches further than htslib does at any version: a header that types
SVLEN as a string rather than as the reserved integer it is. htslib declines
to read such a value into the span at all; the readers here take the number
the writer plainly meant, as they already do for an END the header never
declares. Refusing it would put the pysam-backed reader at odds with the text
ones, which parse the raw field and cannot see a declared type in the first
place.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

#: Symbolic ALT types whose SVLEN counts reference bases, matching htslib's
#: `svlen_on_ref_for_vcf_alt`. An insertion is absent because its SVLEN counts
#: inserted bases rather than reference ones; reading it as a span was an
#: htslib defect, https://github.com/samtools/htslib/issues/1940, and
#: restricting the set is how 1.23 fixed it. That release's own note names
#: only DEL, DUP and CNV: the code includes INV, so a reader following the
#: changelog rather than the source would drop every inversion's span.
_SPANNING_ALT_TYPES = frozenset(("CNV", "DEL", "DUP", "INV"))

#: The range in which htslib keeps an INFO integer, its `BCF_MIN_BT_INT32`
#: and `BCF_MAX_BT_INT32`. Outside it the value is replaced by the missing
#: marker, with the warning "Extreme INFO/SVLEN value encountered and set to
#: missing". The range is not symmetric: the bottom stops eight short of
#: INT32_MIN, a margin reserving BCF's own missing and vector-end markers, so
#: testing a magnitude instead would honour a handful of large negative SVLEN
#: values that htslib discards. Matching it also keeps a colossal value from
#: overflowing a coordinate column rather than merely disagreeing with htslib.
_MIN_SVLEN = -2147483640  # INT32_MIN + 8
_MAX_SVLEN = 2147483647  # INT32_MAX


def _is_reference_spanning_alt(alt: str) -> bool:
    """Whether an SVLEN on this ALT allele measures reference bases."""
    # A sub-type is written `<DEL:ME>`, so the character after the type is
    # either the closing bracket or the separator -- which is what keeps a
    # hypothetical `<DELETION>` out.
    return (
        len(alt) >= 5
        and alt.startswith("<")
        and alt.endswith(">")
        and alt[4] in ">:"
        and alt[1:4] in _SPANNING_ALT_TYPES
    )


def svlen_span(alts: Iterable[str], svlens: Iterable[object]) -> int:
    """Reference bases the largest applicable SVLEN covers, or 0 for none.

    SVLEN carries one value per ALT allele, so values pair with alleles by
    position, and a value against an allele that replaces no reference is not
    a span. The sign is not read: a deletion's length is written negative by
    writers following VCF 4.2 and positive by those following 4.4, and htslib
    takes either as a magnitude.

    The result is one more than that magnitude, because a symbolic allele is
    anchored on the base at POS and describes the bases after it. A value of
    zero yields zero rather than a one-base span, so that a record spelling
    `SVLEN=0` falls back to whatever END and the reference allele say, as it
    does in htslib, and so does one too large for htslib to have believed.
    """
    longest = 0
    # `strict=False`: a writer may give fewer SVLEN values than ALT alleles,
    # or more, and htslib pairs whatever both have rather than refusing.
    for alt, svlen in zip(alts, svlens, strict=False):
        # A declared Integer arrives as an int and an undeclared one as text;
        # anything else -- the None of a '.' value above all -- states no
        # length, and the reference allele still says how far the record
        # reaches.
        if not isinstance(svlen, (int, str)) or not _is_reference_spanning_alt(alt):
            continue
        try:
            value = int(svlen)
        except ValueError:
            continue
        # Range-check the value as written, then take its magnitude: the two
        # do not commute at the bottom of the range.
        if _MIN_SVLEN <= value <= _MAX_SVLEN:
            longest = max(longest, abs(value))
    return longest + 1 if longest else 0
