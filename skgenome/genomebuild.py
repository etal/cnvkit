"""Reference assembly metadata used for sex-chromosome / PAR handling.

The ``GenomeBuild`` value object collects the build-specific data CNVkit needs
to interpret a sample's chromosomes correctly: today, that's the
pseudoautosomal-region (PAR) coordinates on chrX and chrY. Adding a new
supported assembly is a single ``GenomeBuild(...)`` instance in this module.

For callers that already pass a genome name as a string (e.g.
``diploid_parx_genome="grch38"``), :func:`get_genome_build` performs the
lookup and produces a clear ``KeyError`` for unknown names.
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
    """

    name: str
    par_regions: Mapping[str, tuple[int, int]]

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
)

GRCH38 = GenomeBuild(
    name="grch38",
    par_regions={
        "PAR1X": (10000, 2781479),
        "PAR2X": (155701382, 156030895),
        "PAR1Y": (10000, 2781479),
        "PAR2Y": (56887902, 57217415),
    },
)


REGISTERED_BUILDS: dict[str, GenomeBuild] = {
    b.name: b for b in (GRCH37, GRCH38)
}


def get_genome_build(name: str) -> GenomeBuild:
    """Look up a registered genome build by name (case-insensitive).

    Raises ``KeyError`` with a list of supported names if *name* is unknown.
    """
    key = name.lower()
    try:
        return REGISTERED_BUILDS[key]
    except KeyError:
        raise KeyError(
            f"Unknown genome build {name!r}; "
            f"supported: {sorted(REGISTERED_BUILDS)}"
        ) from None


def is_supported_build(name: str) -> bool:
    """Return True if *name* is a registered genome build."""
    return name.lower() in REGISTERED_BUILDS
