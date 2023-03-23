#!/usr/bin/env python
"""Unit tests for the CNVkit library, cnvlib."""
import logging
import unittest

from cnvlib.cnary import CopyNumArray

logging.basicConfig(level=logging.ERROR, format="%(message)s")

# unittest/pomegranate 0.10.0: ImportWarning: can't resolve package from
# __spec__ or __package__, falling back on __name__ and __path__
import warnings

warnings.filterwarnings("ignore", category=ImportWarning)

import numpy as np
from skgenome import GenomicArray, tabio

import cnvlib
from cnvlib import fix, params


def read(*args, **kwargs) -> CopyNumArray:
    return cnvlib.read(*args, **kwargs)


class CNATests(unittest.TestCase):
    """Tests for the CopyNumArray class."""

    def test_empty(self):
        """Instantiate from an empty file."""
        cnarr = cnvlib.read("formats/empty", None)
        self.assertEqual(len(cnarr), 0)

    def test_set_genome_build(self):
        ex_cnr = read("formats/reference-tr.cnn", None)
        self.assertEqual(ex_cnr.genome_build, None, "By default, the genome build is None.")
        self.assertRaises(Exception, ex_cnr.set_genome_build, "doesnt-exist")
        ex_cnr.set_genome_build("grch38")
        self.assertEqual(ex_cnr.genome_build, "grch38", "The genome build can be set.")
        ex_cnr.set_genome_build("GRCh38")
        self.assertEqual(ex_cnr.genome_build, "grch38", "The genome build is converted to lower case.")

    def test_par_and_chrx_filter(self):
        ex_cnr = read("formats/par.cnr", None)
        self.assertFalse(
            ex_cnr.is_par_on_chrx_treated_as_autosomal,
            "By default, chromosome X contains PAR1/2.",
        )
        true_len_literal_x = 10967  # all regions on X
        literal_x = ex_cnr.chromosome == "X"
        self.assertEqual(sum(literal_x), true_len_literal_x)
        filter_x = ex_cnr.chrx_filter()
        self.assertEqual(
            sum(filter_x),
            true_len_literal_x,
            "By default, the filter on chr X returns all of X.",
        )
        self.assertRaises(AssertionError, ex_cnr.par_on_chrx_filter)  # By default, the PAR filter cannot be invoked.

        ex_cnr.treat_par_on_chrx_as_autosomal("grch37")
        par_on_x = ex_cnr.par_on_chrx_filter()
        self.assertEqual(sum(par_on_x), 201, "Once the genome is set, PAR1/2 can be filtered for.")

        par1_overlapping_region_index = 205
        par1_overlapping_region = ex_cnr[par1_overlapping_region_index]
        self.assertEqual(par1_overlapping_region.gene, "PAR1Overlap")
        self.assertTrue(
            par1_overlapping_region.start < params.PAR1_X_GRCH37_END < par1_overlapping_region.end,
            "The region overlaps the PAR1 boundary.",
        )
        self.assertFalse(par_on_x[par1_overlapping_region_index], "The overlapping region is not part of the filter.")

        par2_overlapping_region_index = 10970
        par2_overlapping_region = ex_cnr[par2_overlapping_region_index]
        self.assertEqual(par2_overlapping_region.gene, "PAR2Overlap")
        self.assertTrue(
            par2_overlapping_region.start < params.PAR2_X_GRCH37_START < par2_overlapping_region.end,
            "The region overlaps the PAR2 boundary.",
            )
        self.assertFalse(par_on_x[par2_overlapping_region_index], "The overlapping region is not part of the filter.")

    def test_autosomes(self):
        """Test selection of autosomes specific to CNA. """

        ex_cnr = read("formats/reference-tr.cnn", None)
        auto = ex_cnr.autosomes()
        some_x = (ex_cnr.chromosome == "chrX") & (ex_cnr.end <= 434918)
        some_x_len = some_x.sum()
        self.assertEqual(some_x_len, 3)
        auto_and_some_x = ex_cnr.autosomes(also=some_x)
        self.assertEqual(len(auto_and_some_x), len(auto) + some_x_len)

        ex_cnr.treat_par_on_chrx_as_autosomal(genome_build="grch38")
        auto_with_par_on_chrx = ex_cnr.autosomes()
        len_par = ex_cnr.par_on_chrx_filter().sum()
        self.assertEqual(len_par, 25)
        self.assertEqual(len(auto_with_par_on_chrx), len(auto) + len_par)

        ex_cnr.do_not_treat_par_on_chrx_as_autosomal()
        new_auto = ex_cnr.autosomes()
        self.assertNotEqual(
            len(new_auto), len(auto_with_par_on_chrx), "The CNA can be reset and PAR is within X again."
        )

        # todo: add test par + other also

    def test_basic(self):
        """Test basic container functionality and magic methods."""
        cna = cnvlib.read("formats/reference-tr.cnn", None)
        # Length
        self.assertEqual(len(cna), linecount("formats/reference-tr.cnn") - 1)
        # Equality
        same = cnvlib.read("formats/reference-tr.cnn", None)
        self.assertEqual(cna, same)
        # Item access
        orig = cna[0]
        cna[0] = orig
        cna[3:4] = cna[3:4]
        cna[6:10] = cna[6:10]
        self.assertEqual(tuple(cna[0]), tuple(same[0]))
        self.assertEqual(cna[3:6], same[3:6])

    def test_center_all(self):
        """Test recentering."""
        cna = cnvlib.read("formats/reference-tr.cnn", None)
        # Median-centering an already median-centered array -> no change
        chr1 = cna.in_range("chr1")
        self.assertAlmostEqual(0, np.median(chr1["log2"]), places=1)
        chr1.center_all()
        orig_chr1_cvg = np.median(chr1["log2"])
        self.assertAlmostEqual(0, orig_chr1_cvg)
        # Median-centering resets a shift away from the median
        chr1plus2 = chr1.copy()
        chr1plus2["log2"] += 2.0
        chr1plus2.center_all()
        self.assertAlmostEqual(np.median(chr1plus2["log2"]), orig_chr1_cvg)
        # Other methods for centering are similar for a CN-neutral chromosome
        for method in ("mean", "mode", "biweight"):
            cp = chr1.copy()
            cp.center_all(method)
            self.assertLess(abs(cp["log2"].median() - orig_chr1_cvg), 0.1)

    def test_drop_extra_columns(self):
        """Test removal of optional 'gc' column."""
        cna = cnvlib.read("formats/reference-tr.cnn", None)
        self.assertIn("gc", cna)
        cleaned = cna.drop_extra_columns()
        self.assertNotIn("gc", cleaned)
        self.assertTrue((cleaned["log2"] == cna["log2"]).all())

    def test_guess_xx(self):
        """Guess chromosomal sex from chrX log2 ratio value."""
        for (fname, sample_is_f, ref_is_m) in (
            ("formats/f-on-f.cns", True, False),
            ("formats/f-on-m.cns", True, True),
            ("formats/m-on-f.cns", False, False),
            ("formats/m-on-m.cns", False, True),
            ("formats/amplicon.cnr", False, True),
            ("formats/cl_seq.cns", True, True),
            ("formats/tr95t.cns", True, True),
            ("formats/reference-tr.cnn", False, False),
        ):
            guess = cnvlib.read(fname, None).guess_xx(ref_is_m)
            self.assertEqual(
                guess,
                sample_is_f,
                f"{fname}: guessed XX {guess} but is {sample_is_f}",
            )

    def test_residuals(self):
        cnarr = cnvlib.read("formats/amplicon.cnr", None)
        segments = cnvlib.read("formats/amplicon.cns", None)
        regions = GenomicArray(segments.data).drop_extra_columns()
        for grouping_arg in (None, segments, regions):
            resid = cnarr.residuals(grouping_arg)
            self.assertAlmostEqual(0, resid.mean(), delta=0.3)
            self.assertAlmostEqual(1, np.percentile(resid, 80), delta=0.2)
            self.assertAlmostEqual(2, resid.std(), delta=0.5)

    def test_smooth_log2(self):
        for fname in [
            "formats/amplicon.cnr",
            "formats/wgs-chr17.cnr",
            "formats/p2-20_1.cnr",
            "formats/p2-20_2.cnr",
        ]:
            cnarr = cnvlib.read(fname, None)
            orig_vals = cnarr.log2.values.copy()
            signal = cnarr.smooth_log2()
            self.assertTrue((orig_vals == cnarr.log2.values).all())
            self.assertGreaterEqual(signal.min(), cnarr.log2.min())
            self.assertLessEqual(signal.max(), cnarr.log2.max())


class OtherTests(unittest.TestCase):
    """Tests for other functionality."""

    def test_fix_edge(self):
        """Test the 'edge' bias correction calculations."""
        # NB: With no gap, gain and loss should balance out
        # Wide target, no secondary corrections triggered
        insert_size = 250
        gap_size = np.zeros(1)  # Adjacent
        target_size = np.array([600])
        loss = fix.edge_losses(target_size, insert_size)
        gain = fix.edge_gains(target_size, gap_size, insert_size)
        gain *= 2  # Same on the other side
        self.assertAlmostEqual(loss, gain)
        # Trigger 'loss' correction (target_size < 2 * insert_size)
        target_size = np.array([450])
        self.assertAlmostEqual(
            fix.edge_losses(target_size, insert_size),
            2 * fix.edge_gains(target_size, gap_size, insert_size),
        )
        # Trigger 'gain' correction (target_size + gap_size < insert_size)
        target_size = np.array([300])
        self.assertAlmostEqual(
            fix.edge_losses(target_size, insert_size),
            2 * fix.edge_gains(target_size, gap_size, insert_size),
        )

    # call
    # Test: convert_clonal(x, 1, 2) == convert_diploid(x)


class VATests(unittest.TestCase):
    """Tests for the VariantArray class."""

    def test_read(self):
        """Instantiate from a VCF file."""
        variants = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        self.assertGreater(len(variants), 0)


# == helpers ==


def linecount(filename):
    i = -1
    with open(filename) as handle:
        for i, _line in enumerate(handle):
            pass
        return i + 1


if __name__ == "__main__":
    unittest.main(verbosity=2)
