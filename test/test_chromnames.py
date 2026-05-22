#!/usr/bin/env python
"""Unit tests for skgenome.chromnames."""

import unittest

from skgenome.chromnames import (
    detect_chr_prefix,
    diagnose_missing_chromosome,
    infer_sex_chrom_labels,
    is_alternative_contig,
    is_arabic_numeral_chrom,
    is_autosome,
    is_mitochondrial,
    is_roman_numeral_chrom,
    looks_like_roman_numeral_genome,
)


class IsArabicNumeralChromTests(unittest.TestCase):
    def test_plain_integers(self):
        for name in ("1", "2", "10", "22"):
            self.assertTrue(is_arabic_numeral_chrom(name), name)

    def test_chr_prefix(self):
        for name in ("chr1", "chr22", "Chr3"):
            self.assertTrue(is_arabic_numeral_chrom(name), name)

    def test_sex_and_mito(self):
        for name in ("X", "chrX", "Y", "chrY", "M", "MT", "chrM"):
            self.assertFalse(is_arabic_numeral_chrom(name), name)

    def test_alt_contigs(self):
        for name in ("chr1_random", "Un_gl000220", "HLA-A"):
            self.assertFalse(is_arabic_numeral_chrom(name), name)


class IsRomanNumeralChromTests(unittest.TestCase):
    def test_valid_low(self):
        for name in ("I", "II", "III", "IV", "V", "VI", "IX", "X"):
            self.assertTrue(is_roman_numeral_chrom(name), name)

    def test_valid_yeast(self):
        for name in ("chrI", "chrII", "chrVIII", "chrIX", "chrXI", "chrXVI"):
            self.assertTrue(is_roman_numeral_chrom(name), name)

    def test_invalid_forms(self):
        for name in ("IIII", "VV", "IC", "XM", "VX"):
            self.assertFalse(is_roman_numeral_chrom(name), name)

    def test_not_roman(self):
        for name in ("1", "chr1", "Y", "chrY", "MT", "M", ""):
            self.assertFalse(is_roman_numeral_chrom(name), name)


class IsMitochondrialTests(unittest.TestCase):
    def test_matches(self):
        for name in ("M", "MT", "Mito", "chrM", "chrMT", "ChrMito"):
            self.assertTrue(is_mitochondrial(name), name)

    def test_does_not_match(self):
        for name in ("chrMito_random", "chr1", "X", "MTB", "MTOR"):
            self.assertFalse(is_mitochondrial(name), name)


class IsAlternativeContigTests(unittest.TestCase):
    def test_matches(self):
        for name in (
            "chr1_random",
            "chrUn_gl000220",
            "HLA-A",
            "HLA_DRB1",
            "chr6_alt",
            "chr19_hap1",
            "chr6_hap5",
            "chrEBV",
            "NC_007605",
            "NC-foo",
            "hs37d5_decoy",
        ):
            self.assertTrue(is_alternative_contig(name), name)

    def test_does_not_match(self):
        for name in ("chr1", "chrX", "chrM", "chrI", "chrXVI", "1", "X", "Y", "MT"):
            self.assertFalse(is_alternative_contig(name), name)


class LooksLikeRomanNumeralGenomeTests(unittest.TestCase):
    def test_yeast(self):
        yeast = [f"chr{r}" for r in ("I", "II", "III", "IV", "V", "VIII", "XVI")]
        self.assertTrue(looks_like_roman_numeral_genome(yeast))

    def test_human(self):
        human = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        self.assertFalse(looks_like_roman_numeral_genome(human))

    def test_too_few_roman_names(self):
        # Just chrX in a human dataset is not enough to declare Roman-numeral
        self.assertFalse(looks_like_roman_numeral_genome(["chr1", "chr2", "chrX"]))

    def test_only_ambiguous_single_chars(self):
        # I, V, X alone aren't enough to be sure
        self.assertFalse(looks_like_roman_numeral_genome(["chrI", "chrV", "chrX"]))


class InferSexChromLabelsTests(unittest.TestCase):
    def test_human_chr_prefix(self):
        x, y = infer_sex_chrom_labels(["chr1", "chr2", "chrX", "chrY"])
        self.assertEqual((x, y), ("chrX", "chrY"))

    def test_human_no_prefix(self):
        x, y = infer_sex_chrom_labels(["1", "2", "X", "Y"])
        self.assertEqual((x, y), ("X", "Y"))

    def test_yeast_returns_none(self):
        yeast = [f"chr{r}" for r in ("I", "II", "VIII", "X", "XI", "XVI", "M")]
        self.assertEqual(infer_sex_chrom_labels(yeast), (None, None))

    def test_missing_chromosomes(self):
        # X present, Y absent
        x, y = infer_sex_chrom_labels(["chr1", "chrX"])
        self.assertEqual((x, y), ("chrX", None))

    def test_empty(self):
        self.assertEqual(infer_sex_chrom_labels([]), (None, None))


class DetectChrPrefixTests(unittest.TestCase):
    def test_chr_prefix(self):
        self.assertEqual(detect_chr_prefix(["chr1", "chr2", "chrX"]), "chr")

    def test_no_prefix(self):
        self.assertEqual(detect_chr_prefix(["1", "2", "X"]), "")

    def test_mixed_majority_prefixed(self):
        self.assertEqual(detect_chr_prefix(["chr1", "chr2", "MT"]), "chr")

    def test_empty(self):
        self.assertEqual(detect_chr_prefix([]), "")


class IsAutosomeTests(unittest.TestCase):
    def test_arabic_autosome(self):
        self.assertTrue(is_autosome("chr1", "chrX", "chrY"))
        self.assertTrue(is_autosome("22"))

    def test_sex_excluded(self):
        self.assertFalse(is_autosome("chrX", sex_x_label="chrX", sex_y_label="chrY"))
        self.assertFalse(is_autosome("chrY", sex_x_label="chrX", sex_y_label="chrY"))

    def test_yeast_chr_x_is_autosome(self):
        # In yeast, chrX is autosome 10 -- and infer_sex_chrom_labels returns None
        self.assertTrue(is_autosome("chrX", sex_x_label=None, sex_y_label=None))
        self.assertTrue(is_autosome("chrXVI"))

    def test_mito_is_not_autosome(self):
        # Mito doesn't match either regex (arabic or roman)
        self.assertFalse(is_autosome("chrM"))

    def test_empty(self):
        self.assertFalse(is_autosome(""))


class DiagnoseMissingChromosomeTests(unittest.TestCase):
    def test_suggests_chr_prefix(self):
        msg = diagnose_missing_chromosome("8", ["chr1", "chr2", "chr8", "chrX"])
        self.assertIn("'chr8'", msg)
        self.assertIn("Did you mean", msg)

    def test_suggests_strip_prefix(self):
        msg = diagnose_missing_chromosome("chr8", ["1", "2", "8", "X"])
        self.assertIn("Did you mean '8'", msg)

    def test_present_but_empty(self):
        msg = diagnose_missing_chromosome("chr8", ["chr1", "chr8"])
        self.assertIn("matched 0 rows", msg)
        self.assertIn("NaN", msg)

    def test_no_chromosomes(self):
        msg = diagnose_missing_chromosome("chr1", [])
        self.assertIn("no rows", msg)


if __name__ == "__main__":
    unittest.main()
