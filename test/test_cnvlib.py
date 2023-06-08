#!/usr/bin/env python
"""Unit tests for the CNVkit library, cnvlib."""
import logging
import unittest

logging.basicConfig(level=logging.ERROR, format="%(message)s")

# unittest/pomegranate 0.10.0: ImportWarning: can't resolve package from
# __spec__ or __package__, falling back on __name__ and __path__
import warnings

warnings.filterwarnings("ignore", category=ImportWarning)

import numpy as np
from skgenome import GenomicArray, tabio

import cnvlib
from cnvlib import cnary, fix, params


class CNATests(unittest.TestCase):
    """Tests for the CopyNumArray class."""

    def test_empty(self):
        """Instantiate from an empty file."""
        cnarr = cnvlib.read("formats/empty")
        self.assertEqual(len(cnarr), 0)

    def test_par_and_chrxy_filter(self):
        log2_value = 0.0
        ex_cnr = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 10467178, 10467348, "DFFA", log2_value],
                ["chrX", 640975, 641119, "SHOX", log2_value],  # PAR1X
                ["chrX", 2600000, 2800000, "PAR1X-overlap", log2_value],
                ["chrX", 65727808, 65727914, "MSN", log2_value],
                ["chrX", 155600000, 1557000000, "PAR2X-overlap", log2_value],
                ["chrX", 155767712, 155768156, "SPRY3", log2_value],  # PAR2X
                ["chrY", 640975, 641119, "SHOXY", log2_value],  # PAR1Y
                ["chrY", 2600000, 2800000, "PAR1Y-overlap", log2_value],
                ["chrY", 17877831, 17880784,  "CDY2B", log2_value],
                ["chrY", 56700000, 56900000, "PAR2Y-overlap", log2_value],
                ["chrY", 56950000,56960000, "SPRY3", log2_value],  # PAR2Y
            ]
        )
        # X
        true_number_of_regions_on_x = 5  # all regions on X
        literal_x = ex_cnr[ex_cnr.chromosome == "chrX"]
        self.assertEqual(len(literal_x), true_number_of_regions_on_x, "Baseline assumption is correct.")
        filter_x = ex_cnr[ex_cnr.chr_x_filter()]
        self.assertEqual(len(filter_x), true_number_of_regions_on_x, "By default, the filter on chr X returns all of X.")
        filter_x_diploid_par = ex_cnr[ex_cnr.chr_x_filter("grch38")]
        self.assertTrue(len(filter_x_diploid_par) < true_number_of_regions_on_x, "PAR1/2 can be ignored.")

        par_on_x = ex_cnr[ex_cnr.parx_filter("grch38")]
        self.assertEqual(len(par_on_x), 2, "Filtering for PAR1/2 works fine.")

        par1x_overlapping_region = ex_cnr[2]
        par1x_overlapping_gene = "PAR1X-overlap"
        self.assertEqual(par1x_overlapping_region.gene, par1x_overlapping_gene)
        grch38_par1x_end = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR1X"][1]
        self.assertTrue(par1x_overlapping_region.start < grch38_par1x_end < par1x_overlapping_region.end, "The region overlaps the PAR1 boundary.")
        self.assertTrue(par1x_overlapping_gene not in par_on_x["gene"].tolist(), "The overlapping region is not part of the filter.")

        par2x_overlapping_region = ex_cnr[4]
        par2x_overlapping_gene = "PAR2X-overlap"
        self.assertEqual(par2x_overlapping_region.gene, par2x_overlapping_gene)
        grch38_par2x_start = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR2X"][0]
        self.assertTrue(par2x_overlapping_region.start < grch38_par2x_start < par2x_overlapping_region.end, "The region overlaps the PAR2 boundary.")
        self.assertTrue(par2x_overlapping_gene not in par_on_x["gene"].tolist(), "The overlapping region is not part of the filter.")

        # Y
        true_number_of_regions_on_y = 5  # all regions on y
        literal_y = ex_cnr[ex_cnr.chromosome == "chrY"]
        self.assertEqual(len(literal_y), true_number_of_regions_on_y, "Baseline assumption is correct.")
        filter_y = ex_cnr[ex_cnr.chr_y_filter()]
        self.assertEqual(len(filter_y), true_number_of_regions_on_y, "By default, the filter on chr y returns all of y.")
        filter_y_diploid_par = ex_cnr[ex_cnr.chr_y_filter("grch38")]
        self.assertTrue(len(filter_y_diploid_par) < true_number_of_regions_on_y, "PAR1/2 can be ignored.")

        par_on_y = ex_cnr[ex_cnr.pary_filter("grch38")]
        self.assertEqual(len(par_on_y), 2, "Filtering for PAR1/2 works fine.")

        par1y_overlapping_region = ex_cnr[7]
        par1y_overlapping_gene = "PAR1Y-overlap"
        self.assertEqual(par1y_overlapping_region.gene, par1y_overlapping_gene)
        grch38_par1y_end = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR1Y"][1]
        self.assertTrue(par1y_overlapping_region.start < grch38_par1y_end < par1y_overlapping_region.end, "The region overlaps the PAR1 boundary.")
        self.assertTrue(par1y_overlapping_gene not in par_on_y["gene"].tolist(), "The overlapping region is not part of the filter.")

        par2y_overlapping_region = ex_cnr[9]
        par2y_overlapping_gene = "PAR2Y-overlap"
        self.assertEqual(par2y_overlapping_region.gene, par2y_overlapping_gene)
        grch38_par2y_start = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR2Y"][0]
        self.assertTrue(par2y_overlapping_region.start < grch38_par2y_start < par2y_overlapping_region.end, "The region overlaps the PAR2 boundary.")
        self.assertTrue(par2y_overlapping_gene not in par_on_y["gene"].tolist(), "The overlapping region is not part of the filter.")

    def test_autosomes(self):
        """Test selection of autosomes specific to CNA. """

        ex_cnn = cnvlib.read("formats/par-reference.grch38.cnn")
        auto = ex_cnn.autosomes()
        true_number_of_non_x_non_y_regions = 11
        self.assertEqual(len(auto), true_number_of_non_x_non_y_regions, "Extraction of autosomes works as expected.")

        true_number_of_some_x_regions = 11
        some_x_filter = (ex_cnn.chromosome == "chrX") & (ex_cnn.end <= 154563000)
        some_x = ex_cnn[some_x_filter]
        self.assertEqual(len(some_x), true_number_of_some_x_regions)
        auto_and_some_x = ex_cnn.autosomes(also=some_x_filter)
        self.assertEqual(len(auto_and_some_x), len(auto) + len(some_x), "It is possible to provide further pd.Series as filter.")

        auto_with_parx = ex_cnn.autosomes("grch38")
        parx = ex_cnn[ex_cnn.parx_filter("grch38")]
        self.assertEqual(len(auto_with_parx), true_number_of_non_x_non_y_regions + len(parx), "PAR1/2 is included upon request.")

    def test_basic(self):
        """Test basic container functionality and magic methods."""
        cna = cnvlib.read("formats/reference-tr.cnn")
        # Length
        self.assertEqual(len(cna), linecount("formats/reference-tr.cnn") - 1)
        # Equality
        same = cnvlib.read("formats/reference-tr.cnn")
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
        cna = cnvlib.read("formats/reference-tr.cnn")
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

        # PAR setting influences centering.
        cna1 = cnvlib.read("formats/par-reference.grch38.cnn")
        before1 = np.median(cna1["log2"])
        cna1.center_all()
        after1 = np.median(cna1["log2"])
        cna2 = cnvlib.read("formats/par-reference.grch38.cnn")
        before2 = np.median(cna2["log2"])
        cna2.center_all(diploid_parx_genome="grch38")
        after2 = np.median(cna2["log2"])
        self.assertEqual(before1, before2)
        self.assertNotEqual(after1, after2)

    def test_drop_extra_columns(self):
        """Test removal of optional 'gc' column."""
        cna = cnvlib.read("formats/reference-tr.cnn")
        self.assertIn("gc", cna)
        cleaned = cna.drop_extra_columns()
        self.assertNotIn("gc", cleaned)
        self.assertTrue((cleaned["log2"] == cna["log2"]).all())

    def test_guess_xx(self):
        # TODO: Consider adding a test with --diploid-parx-genome option
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
            guess = cnvlib.read(fname).guess_xx(ref_is_m)
            self.assertEqual(
                guess,
                sample_is_f,
                f"{fname}: guessed XX {guess} but is {sample_is_f}",
            )

    def test_residuals(self):
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segments = cnvlib.read("formats/amplicon.cns")
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
            cnarr = cnvlib.read(fname)
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
