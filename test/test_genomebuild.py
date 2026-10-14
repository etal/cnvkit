#!/usr/bin/env python
"""Unit tests for skgenome.genomebuild."""

import unittest
from dataclasses import FrozenInstanceError

from cnvlib import params
from skgenome.genomebuild import (
    GRCH37,
    GRCH38,
    REGISTERED_BUILDS,
    GenomeBuild,
    get_genome_build,
    is_supported_build,
)


class GenomeBuildTests(unittest.TestCase):
    def test_registered_builds_present(self):
        self.assertIn("grch37", REGISTERED_BUILDS)
        self.assertIn("grch38", REGISTERED_BUILDS)

    def test_par_regions_complete(self):
        for build in (GRCH37, GRCH38):
            for key in ("PAR1X", "PAR2X", "PAR1Y", "PAR2Y"):
                self.assertIn(key, build.par_regions, f"{build.name} missing {key}")
                start, end = build.par_regions[key]
                self.assertLess(start, end, f"{build.name}/{key} bad coords")

    def test_get_genome_build_case_insensitive(self):
        self.assertIs(get_genome_build("GRCh38"), GRCH38)
        self.assertIs(get_genome_build("grch38"), GRCH38)
        self.assertIs(get_genome_build("GRCH38"), GRCH38)

    def test_ucsc_aliases_resolve_to_canonical_build(self):
        """UCSC hg19/hg38 share GRCh37/GRCh38 chrX and chrY coordinates."""
        self.assertIs(get_genome_build("hg19"), GRCH37)
        self.assertIs(get_genome_build("hg38"), GRCH38)
        self.assertIs(get_genome_build("HG38"), GRCH38)
        self.assertEqual(get_genome_build("hg38").name, "grch38")

    def test_get_genome_build_unknown_raises(self):
        # hg18 is NCBI36, whose PAR coordinates differ; it must not resolve.
        with self.assertRaises(KeyError) as cm:
            get_genome_build("hg18")
        # Error message lists accepted names, aliases included, so the user can
        # self-correct.
        self.assertIn("supported", str(cm.exception))
        self.assertIn("hg38", str(cm.exception))

    def test_is_supported_build(self):
        self.assertTrue(is_supported_build("grch37"))
        self.assertTrue(is_supported_build("GRCh38"))
        self.assertTrue(is_supported_build("HG19"))
        self.assertFalse(is_supported_build("hg18"))

    def test_alias_table_is_unambiguous(self):
        """Every name is lowercase and resolves to the build that declares it.

        A collision between any two names, canonical or alias, would silently
        resolve to whichever build came last and apply the wrong PAR
        coordinates; lookup lowercases its argument, so an uppercase alias
        could never match.
        """
        for build in REGISTERED_BUILDS.values():
            for name in (build.name, *build.aliases):
                self.assertEqual(name, name.lower())
                self.assertIs(get_genome_build(name), build)

    def test_frozen(self):
        with self.assertRaises(FrozenInstanceError):
            GRCH37.name = "other"  # type: ignore[misc]

    def test_par_regions_singleton_protected_from_mutation(self):
        """In-place mutation of GenomeBuild.par_regions must be rejected."""
        with self.assertRaises(TypeError):
            GRCH38.par_regions["PAR1X"] = (0, 0)  # type: ignore[index]

    def test_back_compat_via_params(self):
        """cnvlib.params.PSEUDO_AUTSOMAL_REGIONS preserves the historical API."""
        # Every registered build appears as a top-level key, and aliases do not:
        # the historical keys were exactly the canonical names.
        self.assertEqual(set(params.PSEUDO_AUTSOMAL_REGIONS), set(REGISTERED_BUILDS))
        self.assertNotIn("hg38", params.PSEUDO_AUTSOMAL_REGIONS)
        # All four PAR keys round-trip equally for every build
        for build in (GRCH37, GRCH38):
            for key in ("PAR1X", "PAR2X", "PAR1Y", "PAR2Y"):
                self.assertEqual(
                    list(build.par_regions[key]),
                    params.PSEUDO_AUTSOMAL_REGIONS[build.name][key],
                )
        # Re-exported inner values are mutable lists, not tuples or proxies
        # (so callers that mutated the legacy list-of-ints keep working).
        coords = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR1X"]
        self.assertIsInstance(coords, list)
        # And the supported-genomes view stays consistent
        self.assertEqual(
            set(params.SUPPORTED_GENOMES_FOR_PAR_HANDLING),
            set(params.PSEUDO_AUTSOMAL_REGIONS),
        )


if __name__ == "__main__":
    unittest.main()
