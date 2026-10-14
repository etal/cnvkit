"""Reference assembly metadata used for sex-chromosome / PAR handling.

The ``GenomeBuild`` value object collects the build-specific data CNVkit needs
to interpret a sample's chromosomes correctly: today, that's the
pseudoautosomal-region (PAR) coordinates on chrX and chrY. Adding a new
supported assembly is a single ``GenomeBuild(...)`` instance in this module.

For callers that already pass a genome name as a string (e.g.
``diploid_parx_genome="grch38"``), :func:`get_genome_build` performs the
lookup, by canonical name or alias, and produces a clear ``KeyError`` for
unknown names.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping


@dataclass(frozen=True)
class GenomeBuild:
    """Immutable description of a reference assembly's PAR coordinates.

    ``par_regions`` is keyed by ``"PAR1X"``, ``"PAR2X"``, ``"PAR1Y"``,
    ``"PAR2Y"`` and gives the half-open coordinate range as a tuple. The
    backing mapping is wrapped in ``MappingProxyType`` on construction so
    that the module-level singletons (``GRCH37``, ``GRCH38``) cannot be
    mutated in place by callers.

    ``aliases`` lists further names :func:`get_genome_build` accepts for this
    build: other names for the assembly whose chromosome X and Y coordinates
    are identical, so the same PAR coordinates apply. UCSC hg19 and hg38 are
    aliases of GRCh37 and GRCh38 in this sense. Aliases must be lowercase,
    since lookup lowercases its argument.
    """

    name: str
    par_regions: Mapping[str, tuple[int, int]]
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        # Defend the singleton against in-place mutation of the backing dict.
        object.__setattr__(
            self, "par_regions", MappingProxyType(dict(self.par_regions))
        )


GRCH37 = GenomeBuild(
    name="grch37",
    par_regions={
        "PAR1X": (60000, 2699520),
        "PAR2X": (154931043, 155260560),
        "PAR1Y": (10000, 2649520),
        "PAR2Y": (59034049, 59363566),
    },
    aliases=("hg19",),
)

GRCH38 = GenomeBuild(
    name="grch38",
    par_regions={
        "PAR1X": (10000, 2781479),
        "PAR2X": (155701382, 156030895),
        "PAR1Y": (10000, 2781479),
        "PAR2Y": (56887902, 57217415),
    },
    aliases=("hg38",),
)


REGISTERED_BUILDS: dict[str, GenomeBuild] = {b.name: b for b in (GRCH37, GRCH38)}

# Every name lookup accepts: canonical names plus aliases. REGISTERED_BUILDS
# keeps one entry per assembly, because cnvlib.params re-exports its keys as
# the historical PSEUDO_AUTSOMAL_REGIONS table.
_BUILDS_BY_NAME: dict[str, GenomeBuild] = {
    name: build
    for build in REGISTERED_BUILDS.values()
    for name in (build.name, *build.aliases)
}


def get_genome_build(name: str) -> GenomeBuild:
    """Look up a registered genome build by name (case-insensitive).

    Raises ``KeyError`` with a list of supported names if *name* is unknown.
    """
    key = name.lower()
    try:
        return _BUILDS_BY_NAME[key]
    except KeyError:
        raise KeyError(
            f"Unknown genome build {name!r}; supported: {sorted(_BUILDS_BY_NAME)}"
        ) from None


def is_supported_build(name: str) -> bool:
    """Return True if *name* is a registered genome build or alias."""
    return name.lower() in _BUILDS_BY_NAME
