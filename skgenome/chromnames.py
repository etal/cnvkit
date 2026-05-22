"""Chromosome name classification and inference.

This module centralizes the heuristics CNVkit uses to interpret chromosome
names across reference assemblies: detecting autosomes vs. sex chromosomes
vs. mitochondria vs. alternative contigs, recognizing the "chr" prefix
convention, and inferring which label a dataset uses for the X and Y
chromosomes.

Design notes
------------

These functions split into two groups:

* Pure single-string predicates (``is_arabic_numeral_chrom``,
  ``is_roman_numeral_chrom``, ``is_mitochondrial``,
  ``is_alternative_contig``) classify a name in isolation.
* Context-aware classifiers (``looks_like_roman_numeral_genome``,
  ``infer_sex_chrom_labels``, ``is_autosome``) consider the full set of
  chromosome names in a dataset, which is necessary to disambiguate cases
  like ``chrX`` (a sex chromosome in mammals, but Roman numeral 10 in
  yeast).

When a dataset's chromosome names look unfamiliar, the public callers
fall back to a permissive interpretation rather than silently discarding
data. See ``GenomicArray.autosomes`` for an example of that policy.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


_ARABIC_NUMERAL = re.compile(r"^(?:chr|Chr)?\d{1,3}$")

_ALTERNATIVE_CONTIG = re.compile(
    "|".join(
        (
            r"^chrEBV$",
            r"^NC[_\-]",
            r"_random$",
            r"Un_",
            r"^HLA[\-_]",
            r"_alt$",
            r"hap\d+$",
            r"_decoy$",
        )
    )
)

_MITOCHONDRIAL = re.compile(r"^(?:chr|Chr)?(?:M|MT|Mito)$", re.IGNORECASE)


def _to_roman(n: int) -> str:
    vals = [
        (100, "C"),
        (90, "XC"),
        (50, "L"),
        (40, "XL"),
        (10, "X"),
        (9, "IX"),
        (5, "V"),
        (4, "IV"),
        (1, "I"),
    ]
    out: list[str] = []
    for v, sym in vals:
        while n >= v:
            out.append(sym)
            n -= v
    return "".join(out)


# 1..49 covers every known Roman-numeral chromosome assembly (yeast has 16,
# fission yeast 8, other fungi up to ~30). Using a frozenset built from valid
# Roman numerals rejects invalid forms like "IIII" or "VV" without a regex.
_VALID_ROMAN_NUMERALS = frozenset(_to_roman(i) for i in range(1, 50))


def _strip_chr_prefix(name: str) -> str:
    """Return *name* with a leading 'chr' or 'Chr' removed, if present."""
    return name.removeprefix("chr").removeprefix("Chr")


def is_arabic_numeral_chrom(name: str) -> bool:
    """Match plain integer chromosome names like ``1`` or ``chr22``."""
    return bool(_ARABIC_NUMERAL.match(name))


def is_roman_numeral_chrom(name: str) -> bool:
    """Match Roman-numeral chromosome names like ``chrI`` or ``XVI``.

    Only valid Roman numerals up to L (50) are accepted; invalid forms
    such as ``IIII`` or ``VV`` return False.
    """
    stripped = _strip_chr_prefix(name)
    return stripped in _VALID_ROMAN_NUMERALS


def is_mitochondrial(name: str) -> bool:
    """Match mitochondrial chromosome names: ``M``, ``MT``, ``Mito``, etc."""
    return bool(_MITOCHONDRIAL.match(name))


def is_alternative_contig(name: str) -> bool:
    """Match alternative or unplaced contigs: ``_random``, ``Un_``, ``HLA-``, etc."""
    return bool(_ALTERNATIVE_CONTIG.search(name))


def looks_like_roman_numeral_genome(names: Iterable[str]) -> bool:
    """Heuristic: does this chromosome set use Roman-numeral autosomes?

    A human dataset with just ``chrX`` and ``chrI`` (a typo or odd subset)
    is not enough evidence — we require at least three valid Roman-numeral
    names AND at least one unambiguous name (more than one Roman character)
    to avoid false positives from single-letter names ``I``/``V``/``X`` that
    can coincide with normal human chromosomes.
    """
    roman_names = [n for n in names if is_roman_numeral_chrom(n)]
    if len(roman_names) < 3:
        return False
    return any(len(_strip_chr_prefix(n)) > 1 for n in roman_names)


def infer_sex_chrom_labels(
    names: Iterable[str],
) -> tuple[str | None, str | None]:
    """Detect the X and Y chromosome labels used in a chromosome set.

    Returns a ``(x_label, y_label)`` pair, each of which may be None.

    Returns ``(None, None)`` when the chromosome set looks like a Roman-numeral
    genome (e.g. yeast), where ``chrX`` is autosome 10 rather than a sex
    chromosome.
    """
    name_set = {str(n) for n in names}
    if not name_set:
        return None, None
    if looks_like_roman_numeral_genome(name_set):
        return None, None
    x_label = next((c for c in ("chrX", "X") if c in name_set), None)
    y_label = next((c for c in ("chrY", "Y") if c in name_set), None)
    return x_label, y_label


def detect_chr_prefix(names: Iterable[str]) -> str:
    """Return ``"chr"`` if most chromosome names use that prefix, else ``""``.

    Inspects up to the first 32 names to keep this cheap on large sets.
    """
    sample = []
    for i, name in enumerate(names):
        if i >= 32:
            break
        sample.append(name)
    if not sample:
        return ""
    n_prefixed = sum(1 for n in sample if n.startswith("chr"))
    return "chr" if n_prefixed * 2 > len(sample) else ""


def is_autosome(
    name: str,
    sex_x_label: str | None = None,
    sex_y_label: str | None = None,
) -> bool:
    """Test whether *name* is a standard autosome.

    Accepts arabic-numeral names (``1``, ``chr1``) and Roman-numeral names
    (``chrI``, ``XVI``). When *sex_x_label* or *sex_y_label* is supplied,
    those names are excluded so that e.g. ``chrX`` is recognized as a sex
    chromosome rather than Roman numeral 10.
    """
    if not name:
        return False
    if sex_x_label is not None and name == sex_x_label:
        return False
    if sex_y_label is not None and name == sex_y_label:
        return False
    return is_arabic_numeral_chrom(name) or is_roman_numeral_chrom(name)


def diagnose_missing_chromosome(
    requested: str,
    available: Iterable[str],
) -> str:
    """Build a user-facing message for "0 probes" / "chromosome absent" cases.

    Suggests a ``chr``/no-``chr`` prefix variant if one would match, and
    summarizes which chromosomes are present.
    """
    available_sorted = sorted({str(c) for c in available})
    if not available_sorted:
        return f"The data contains no rows; cannot show region {requested!r}."
    sample = ", ".join(repr(c) for c in available_sorted[:8])
    if len(available_sorted) > 8:
        sample += f", ... (+{len(available_sorted) - 8} more)"
    if requested in available_sorted:
        return (
            f"Region {requested!r} matched 0 rows; "
            f"check whether log2 values were dropped as NaN."
        )
    alt = (
        requested.removeprefix("chr")
        if requested.startswith("chr")
        else "chr" + requested
    )
    suggestion = f" Did you mean {alt!r}?" if alt in available_sorted else ""
    return (
        f"Chromosome {requested!r} is not present in the data. "
        f"Available chromosomes: {sample}.{suggestion}"
    )
