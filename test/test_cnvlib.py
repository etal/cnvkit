#!/usr/bin/env python
"""Unit tests for the CNVkit library, cnvlib."""

import logging
import math
import unittest
import warnings

import pytest

logging.basicConfig(level=logging.ERROR, format="%(message)s")

import numpy as np
import pandas as pd
from skgenome import GenomicArray, tabio

import cnvlib
from cnvlib import (
    cnary,
    fix,
    params,
    plots,
    reports,
    segfilters,
    segmentation,
    segmetrics,
    vary,
)
from conftest import linecount


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
                ["chrY", 17877831, 17880784, "CDY2B", log2_value],
                ["chrY", 56700000, 56900000, "PAR2Y-overlap", log2_value],
                ["chrY", 56950000, 56960000, "SPRY3", log2_value],  # PAR2Y
            ]
        )
        # X
        true_number_of_regions_on_x = 5  # all regions on X
        literal_x = ex_cnr[ex_cnr.chromosome == "chrX"]
        self.assertEqual(
            len(literal_x),
            true_number_of_regions_on_x,
            "Baseline assumption is correct.",
        )
        filter_x = ex_cnr[ex_cnr.chr_x_filter()]
        self.assertEqual(
            len(filter_x),
            true_number_of_regions_on_x,
            "By default, the filter on chr X returns all of X.",
        )
        filter_x_diploid_par = ex_cnr[ex_cnr.chr_x_filter("grch38")]
        self.assertTrue(
            len(filter_x_diploid_par) < true_number_of_regions_on_x,
            "PAR1/2 can be ignored.",
        )

        par_on_x = ex_cnr[ex_cnr.parx_filter("grch38")]
        self.assertEqual(len(par_on_x), 2, "Filtering for PAR1/2 works fine.")

        par1x_overlapping_region = ex_cnr[2]
        par1x_overlapping_gene = "PAR1X-overlap"
        self.assertEqual(par1x_overlapping_region.gene, par1x_overlapping_gene)
        grch38_par1x_end = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR1X"][1]
        self.assertTrue(
            par1x_overlapping_region.start
            < grch38_par1x_end
            < par1x_overlapping_region.end,
            "The region overlaps the PAR1 boundary.",
        )
        self.assertTrue(
            par1x_overlapping_gene not in par_on_x["gene"].tolist(),
            "The overlapping region is not part of the filter.",
        )

        par2x_overlapping_region = ex_cnr[4]
        par2x_overlapping_gene = "PAR2X-overlap"
        self.assertEqual(par2x_overlapping_region.gene, par2x_overlapping_gene)
        grch38_par2x_start = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR2X"][0]
        self.assertTrue(
            par2x_overlapping_region.start
            < grch38_par2x_start
            < par2x_overlapping_region.end,
            "The region overlaps the PAR2 boundary.",
        )
        self.assertTrue(
            par2x_overlapping_gene not in par_on_x["gene"].tolist(),
            "The overlapping region is not part of the filter.",
        )

        # Y
        true_number_of_regions_on_y = 5  # all regions on y
        literal_y = ex_cnr[ex_cnr.chromosome == "chrY"]
        self.assertEqual(
            len(literal_y),
            true_number_of_regions_on_y,
            "Baseline assumption is correct.",
        )
        filter_y = ex_cnr[ex_cnr.chr_y_filter()]
        self.assertEqual(
            len(filter_y),
            true_number_of_regions_on_y,
            "By default, the filter on chr y returns all of y.",
        )
        filter_y_diploid_par = ex_cnr[ex_cnr.chr_y_filter("grch38")]
        self.assertTrue(
            len(filter_y_diploid_par) < true_number_of_regions_on_y,
            "PAR1/2 can be ignored.",
        )

        par_on_y = ex_cnr[ex_cnr.pary_filter("grch38")]
        self.assertEqual(len(par_on_y), 2, "Filtering for PAR1/2 works fine.")

        par1y_overlapping_region = ex_cnr[7]
        par1y_overlapping_gene = "PAR1Y-overlap"
        self.assertEqual(par1y_overlapping_region.gene, par1y_overlapping_gene)
        grch38_par1y_end = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR1Y"][1]
        self.assertTrue(
            par1y_overlapping_region.start
            < grch38_par1y_end
            < par1y_overlapping_region.end,
            "The region overlaps the PAR1 boundary.",
        )
        self.assertTrue(
            par1y_overlapping_gene not in par_on_y["gene"].tolist(),
            "The overlapping region is not part of the filter.",
        )

        par2y_overlapping_region = ex_cnr[9]
        par2y_overlapping_gene = "PAR2Y-overlap"
        self.assertEqual(par2y_overlapping_region.gene, par2y_overlapping_gene)
        grch38_par2y_start = params.PSEUDO_AUTSOMAL_REGIONS["grch38"]["PAR2Y"][0]
        self.assertTrue(
            par2y_overlapping_region.start
            < grch38_par2y_start
            < par2y_overlapping_region.end,
            "The region overlaps the PAR2 boundary.",
        )
        self.assertTrue(
            par2y_overlapping_gene not in par_on_y["gene"].tolist(),
            "The overlapping region is not part of the filter.",
        )

    def test_gene_coords_by_name(self):
        """`-g` labels only requested genes, not co-binned neighbors (gh#458).

        A bin whose `gene` column packs several names (e.g. "ERBB2,MIR4728")
        must not surface the unrequested neighbor (MIR4728) when only ERBB2
        is selected.
        """
        cnarr = cnary.CopyNumArray.from_rows(
            [
                ["chr17", 37800000, 37850000, "STARD3", 0.0],
                ["chr17", 37850000, 37860000, "ERBB2,MIR4728", 0.0],
                ["chr17", 37860000, 37870000, "ERBB2", 0.0],
                ["chr17", 37880000, 37890000, "GRB7", 0.0],
            ]
        )
        # Single requested gene: co-binned MIR4728 must be hidden
        coords = plots.gene_coords_by_name(cnarr, ["ERBB2"])
        self.assertEqual(list(coords.keys()), ["chr17"])
        self.assertEqual(len(coords["chr17"]), 1)
        start, end, name = coords["chr17"][0]
        self.assertEqual(name, "ERBB2")
        self.assertEqual((start, end), (37850000, 37870000))
        # When MIR4728 *is* requested, it should appear
        names_seen = set()
        for _s, _e, label in plots.gene_coords_by_name(
            cnarr, ["ERBB2", "MIR4728"]
        )["chr17"]:
            names_seen.update(label.split(","))
        self.assertEqual(names_seen, {"ERBB2", "MIR4728"})

    def test_autosomes(self):
        """Test selection of autosomes specific to CNA."""

        ex_cnn = cnvlib.read("formats/par-reference.grch38.cnn")
        auto = ex_cnn.autosomes()
        true_number_of_non_x_non_y_regions = 11
        self.assertEqual(
            len(auto),
            true_number_of_non_x_non_y_regions,
            "Extraction of autosomes works as expected.",
        )

        true_number_of_some_x_regions = 11
        some_x_filter = (ex_cnn.chromosome == "chrX") & (ex_cnn.end <= 154563000)
        some_x = ex_cnn[some_x_filter]
        self.assertEqual(len(some_x), true_number_of_some_x_regions)
        auto_and_some_x = ex_cnn.autosomes(also=some_x_filter)
        self.assertEqual(
            len(auto_and_some_x),
            len(auto) + len(some_x),
            "It is possible to provide further pd.Series as filter.",
        )

        auto_with_parx = ex_cnn.autosomes("grch38")
        parx = ex_cnn[ex_cnn.parx_filter("grch38")]
        self.assertEqual(
            len(auto_with_parx),
            true_number_of_non_x_non_y_regions + len(parx),
            "PAR1/2 is included upon request.",
        )

    def test_yeast_chromosomes(self):
        """Roman-numeral (yeast) chromosomes are handled correctly.

        Regression test for GitHub issue #889.
        """
        log2_value = 0.0
        rows = [
            (f"chr{r}", 100, 200, "gene", log2_value)
            for r in ("I", "II", "VIII", "IX", "X", "XI", "XVI")
        ] + [("chrM", 100, 200, "gene", log2_value)]
        yeast = cnary.CopyNumArray.from_rows(rows)
        # No sex chromosomes detected -- chrX is autosome 10 in yeast.
        self.assertIsNone(yeast.chr_x_label)
        self.assertIsNone(yeast.chr_y_label)
        # Sex-chromosome filters return all-False, not exceptions.
        self.assertFalse(yeast.chr_x_filter().any())
        self.assertFalse(yeast.chr_y_filter().any())
        # autosomes() includes all 7 Roman-numeral chromosomes (incl. chrX)
        # but excludes chrM.
        auto = yeast.autosomes()
        self.assertEqual(len(auto), 7)
        self.assertNotIn("chrM", set(auto.chromosome))
        self.assertIn("chrX", set(auto.chromosome))
        # guess_xx returns None when there's no X/Y to compare.
        self.assertIsNone(yeast.guess_xx(verbose=False))

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
        for fname, sample_is_f, ref_is_m in (
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

    @pytest.mark.slow
    def test_smooth_log2(self):
        for fname in [
            "formats/amplicon.cnr",
            "formats/wgs-chr17.cnr",
            "formats/p2-20_1.cnr",
            "formats/p2-20_2.cnr",
        ]:
            cnarr = cnvlib.read(fname)
            orig_vals = cnarr.log2.to_numpy().copy()
            signal = cnarr.smooth_log2()
            self.assertTrue((orig_vals == cnarr.log2.to_numpy()).all())
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

    def test_center_by_window_single_bin(self):
        """center_by_window survives one surviving bin (#891).

        Near-zero antitarget coverage can leave a single bin after the
        bad-bin/NaN filters in load_adjust_coverages. The bias-correction
        smoothing (rolling_median) then received len-1 input, where the
        fraction ``max(0.01, len ** -0.5) == 1.0`` is an invalid width that
        crashed ``_width2wing`` (and len-0 tripped its ``wing >= 1`` assert),
        terminating ``batch`` silently without ever writing a .cnr.
        """
        one_bin = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"],
                    "start": [1000],
                    "end": [2000],
                    "gene": ["G"],
                    "log2": [0.42],
                    "depth": [50.0],
                }
            )
        )
        # frac as computed in load_adjust_coverages for one surviving bin
        frac = max(0.01, len(one_bin) ** -0.5)
        self.assertEqual(frac, 1.0)
        result = fix.center_by_window(one_bin, frac, np.array([0.5]))
        self.assertEqual(len(result), 1)
        self.assertFalse(np.isnan(result["log2"]).any())

    # call
    # Test: convert_clonal(x, 1, 2) == convert_diploid(x)

    def test_segment_mean_nan_weights(self):
        """segment_mean handles NaN weights correctly."""
        cnarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 4,
                    "start": [0, 1000, 2000, 3000],
                    "end": [1000, 2000, 3000, 4000],
                    "gene": ["G"] * 4,
                    "log2": [0.2, 0.4, 0.6, 0.8],
                    "weight": [1.0, 2.0, 1.0, 2.0],
                }
            )
        )
        # Baseline: weighted average of all bins
        expected_clean = np.average([0.2, 0.4, 0.6, 0.8], weights=[1, 2, 1, 2])
        self.assertAlmostEqual(segmetrics.segment_mean(cnarr), expected_clean)

        # Partial NaN: weighted average of valid-weight bins only
        cnarr_nan = cnarr.copy()
        cnarr_nan["weight"] = [np.nan, 2.0, np.nan, 2.0]
        expected_partial = np.average([0.4, 0.8], weights=[2, 2])
        self.assertAlmostEqual(segmetrics.segment_mean(cnarr_nan), expected_partial)

        # All NaN: falls through to unweighted mean
        cnarr_allnan = cnarr.copy()
        cnarr_allnan["weight"] = np.nan
        expected_allnan = np.mean([0.2, 0.4, 0.6, 0.8])
        self.assertAlmostEqual(segmetrics.segment_mean(cnarr_allnan), expected_allnan)

    def test_apply_weights_no_nan(self):
        """apply_weights never produces NaN weight values."""
        n = 20
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * n,
                "start": np.arange(0, n * 1000, 1000),
                "end": np.arange(1000, n * 1000 + 1000, 1000),
                "gene": ["GeneA"] * n,
                "log2": np.random.default_rng(42).normal(0, 0.3, n),
                "depth": np.random.default_rng(42).uniform(50, 200, n),
            }
        )
        cnarr = cnary.CopyNumArray(data)
        ref_data = data.copy()
        ref_data["spread"] = 0.3
        ref_matched = cnary.CopyNumArray(ref_data)

        # Normal case — weights should all be valid
        result = fix.apply_weights(cnarr, ref_matched, "log2", "spread")
        self.assertFalse(np.isnan(result["weight"]).any())

        # NaN spread — should produce valid weights, not NaN
        ref_nan = ref_matched.copy()
        ref_nan["spread"] = np.nan
        result = fix.apply_weights(cnarr, ref_nan, "log2", "spread")
        self.assertFalse(np.isnan(result["weight"]).any())

    @staticmethod
    def _cna_with_nan_log2():
        """A CopyNumArray whose middle bins carry NaN / sentinel log2 values."""
        log2 = [0.1, np.nan, params.NULL_LOG2_COVERAGE, 0.2, np.nan, -0.1]
        n = len(log2)
        return cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(0, n * 1000, 1000),
                    "end": np.arange(1000, n * 1000 + 1000, 1000),
                    "gene": ["G"] * n,
                    "log2": log2,
                    "depth": [100.0, 100.0, 0.0, 100.0, np.nan, 100.0],
                    "weight": [1.0] * n,
                }
            )
        )

    def test_drop_low_coverage_drops_nan(self):
        """drop_low_coverage removes NaN-log2 bins, not just the sentinel (#521).

        A bare ``log2 < min_cvg`` comparison is False for NaN, so without
        explicit NaN handling, NaN bins survive into segmentation and trigger
        the 'invalid value encountered in greater' RuntimeWarning / CBS crash.
        """
        cnarr = self._cna_with_nan_log2()
        kept = cnarr.drop_low_coverage()
        # No NaN and no sentinel survives
        self.assertFalse(np.isnan(kept["log2"]).any())
        self.assertTrue((kept["log2"] > params.NULL_LOG2_COVERAGE).all())
        # The three good bins remain (0.1, 0.2, -0.1)
        self.assertEqual(len(kept), 3)

    def test_mask_bad_bins_flags_nan(self):
        """mask_bad_bins treats NaN log2/spread reference bins as bad (#521)."""
        n = 5
        ref = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(0, n * 1000, 1000),
                    "end": np.arange(1000, n * 1000 + 1000, 1000),
                    "gene": ["G"] * n,
                    "log2": [0.0, np.nan, 0.0, 0.0, 0.0],
                    "spread": [0.1, 0.1, np.nan, 0.1, 0.1],
                    "depth": [100.0, 100.0, 100.0, 100.0, 100.0],
                }
            )
        )
        mask = fix.mask_bad_bins(ref)
        # NaN-log2 (idx 1) and NaN-spread (idx 2) bins must be flagged bad
        self.assertTrue(mask.iloc[1])
        self.assertTrue(mask.iloc[2])
        # The clean bins are kept
        self.assertFalse(mask.iloc[0])
        self.assertFalse(mask.iloc[3])
        self.assertFalse(mask.iloc[4])

    def test_transfer_fields_nan_gene(self):
        """transfer_fields handles NaN gene names without crashing (issue #900)."""
        n = 10
        cnarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(0, n * 1000, 1000),
                    "end": np.arange(1000, n * 1000 + 1000, 1000),
                    "gene": ["GeneA", float("nan"), "GeneB"] * 3 + [float("nan")],
                    "log2": np.zeros(n),
                    "depth": np.ones(n) * 100.0,
                    "weight": np.ones(n),
                }
            )
        )
        # One segment covering all bins
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"],
                    "start": [0],
                    "end": [n * 1000],
                    "gene": ["-"],
                    "log2": [0.0],
                    "probes": [n],
                    "weight": [0.0],
                }
            )
        )
        result = segmentation.transfer_fields(segarr, cnarr)
        gene_val = result["gene"].iat[0]
        self.assertIsInstance(gene_val, str)
        self.assertNotIn("nan", gene_val.lower())
        self.assertIn("GeneA", gene_val)
        self.assertIn("GeneB", gene_val)

    def test_transfer_fields_genes_near_segment_end(self):
        """Genes from bins near a segment's end are kept in its label (#688).

        Uses the reported EGFR geometry: bins at chr7:55,018,770-55,019,423
        fall inside the segment chr7:54,246,732-55,031,592 (well before its
        end), and more EGFR bins fall in the adjacent segment. EGFR must appear
        in *both* segment labels -- the original report dropped it from the
        first, where the gene-bearing bins sit near the segment's end.
        """
        bins = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr7"] * 6,
                    "start": [54246732, 54800000, 55018770, 55019096,
                              55032092, 55088166],
                    "end": [54300000, 54900000, 55019096, 55019423,
                            55032193, 55088469],
                    "gene": ["VSTM2A", "SEC61G", "EGFR", "EGFR",
                             "EGFR", "EGFR"],
                    "log2": np.zeros(6),
                    "depth": np.ones(6) * 100.0,
                    "weight": np.ones(6),
                }
            )
        )
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr7", "chr7"],
                    "start": [54246732, 55032092],
                    "end": [55031592, 55365525],
                    "gene": ["-", "-"],
                    "log2": [0.0, 0.0],
                    "probes": [4, 2],
                    "weight": [0.0, 0.0],
                }
            )
        )
        result = segmentation.transfer_fields(segarr, bins)
        self.assertIn("EGFR", result["gene"].iat[0])
        self.assertIn("VSTM2A", result["gene"].iat[0])
        self.assertIn("EGFR", result["gene"].iat[1])

    def test_transfer_fields_nan_weights(self):
        """transfer_fields handles NaN bin weights without NaN in .cns output."""
        n = 10
        cnarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(0, n * 1000, 1000),
                    "end": np.arange(1000, n * 1000 + 1000, 1000),
                    "gene": ["GeneA"] * n,
                    "log2": np.zeros(n),
                    "depth": np.ones(n) * 100.0,
                    "weight": [
                        1.0,
                        np.nan,
                        1.0,
                        np.nan,
                        1.0,
                        1.0,
                        np.nan,
                        1.0,
                        1.0,
                        1.0,
                    ],
                }
            )
        )
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"],
                    "start": [0],
                    "end": [n * 1000],
                    "gene": ["-"],
                    "log2": [0.0],
                    "probes": [n],
                    "weight": [0.0],
                }
            )
        )
        result = segmentation.transfer_fields(segarr, cnarr)
        # Segment weight must not be NaN
        self.assertFalse(np.isnan(result["weight"].iat[0]))
        # Should be the sum of non-NaN weights (7.0)
        self.assertAlmostEqual(result["weight"].iat[0], 7.0)
        # Depth should be valid (weighted mean of non-NaN bins)
        self.assertFalse(np.isnan(result["depth"].iat[0]))

        # All-NaN weights: segment weight should be 0, depth should be 0
        cnarr_allnan = cnarr.copy()
        cnarr_allnan["weight"] = np.nan
        result2 = segmentation.transfer_fields(segarr.copy(), cnarr_allnan)
        self.assertEqual(result2["weight"].iat[0], 0.0)
        self.assertEqual(result2["depth"].iat[0], 0.0)

    def test_do_segmentation_drops_nan_log2(self):
        """do_segmentation tolerates NaN-log2 bins on the default path (#881).

        Without --drop-low-coverage (skip_low=False), NaN-log2 bins are never
        filtered before drop_outliers' Savitzky-Golay smoother, whose scipy
        lstsq rejects non-finite input ("array must not contain infs or NaNs"),
        crashing segmentation long before reaching DNAcopy. The NaN bins must be
        dropped first. Uses the pure-Python 'haar' method so no R is required;
        the savgol outlier path that crashed is shared by every method.
        """
        n = 120  # > drop_outliers' width (50) so savgol actually runs
        rng = np.random.default_rng(0)
        log2 = rng.normal(0, 0.2, n)
        log2[5] = np.nan  # inside the leading savgol edge window
        log2[60] = np.nan
        cnarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(0, n * 1000, 1000),
                    "end": np.arange(1000, n * 1000 + 1000, 1000),
                    "gene": ["G"] * n,
                    "log2": log2,
                    "depth": np.ones(n) * 100.0,
                    "weight": np.ones(n),
                }
            )
        )
        # Must not raise, and no NaN should leak into the segment log2 values
        cns = segmentation.do_segmentation(cnarr, "haar", processes=1)
        self.assertGreater(len(cns), 0)
        self.assertFalse(np.isnan(cns["log2"].to_numpy()).any())

    def test_squash_region_nan_weights(self):
        """squash_region handles NaN weights without NaN in output."""
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1", "chr1"],
                "start": [0, 1000, 2000],
                "end": [1000, 2000, 3000],
                "gene": ["G1", "G2", "G3"],
                "log2": [0.1, 0.2, 0.3],
                "probes": [5, 5, 5],
                "weight": [1.0, np.nan, 1.0],
            }
        )
        result = segfilters.squash_region(df)
        self.assertFalse(np.isnan(result["weight"].iat[0]))
        self.assertAlmostEqual(result["weight"].iat[0], 2.0)
        # log2 should be weighted average of bins with valid weights
        self.assertAlmostEqual(
            result["log2"].iat[0], np.average([0.1, 0.3], weights=[1.0, 1.0])
        )

        # All-NaN weights: should fall through to unweighted mean
        df_allnan = df.copy()
        df_allnan["weight"] = np.nan
        result2 = segfilters.squash_region(df_allnan)
        self.assertAlmostEqual(result2["log2"].iat[0], np.mean([0.1, 0.2, 0.3]))

    def test_bic_nan_weights(self):
        """BIC filter handles NaN segment weights without NaN RSS."""
        df = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [0, 1000, 2000, 3000, 4000],
                "end": [1000, 2000, 3000, 4000, 5000],
                "gene": ["A", "B", "C", "D", "E"],
                "log2": [0.0, 0.1, 0.0, -0.1, 0.0],
                "probes": [10, 10, 10, 10, 10],
                "weight": [1.0, np.nan, 2.0, np.nan, 1.0],
            }
        )
        segarr = cnary.CopyNumArray(df)
        # Should not crash; NaN-weight segments get fallback RSS
        result = segfilters.bic(segarr)
        self.assertGreater(len(result), 0)
        self.assertFalse(np.isnan(result["weight"]).any())

    def test_compute_gene_stats_nan_weights(self):
        """compute_gene_stats handles NaN weights for CI/PI without crashing."""
        bins = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 5,
                    "start": np.arange(0, 5000, 1000),
                    "end": np.arange(1000, 6000, 1000),
                    "gene": ["G"] * 5,
                    "log2": [0.1, 0.2, 0.1, 0.3, 0.2],
                    "weight": [1.0, np.nan, 1.0, np.nan, 1.0],
                }
            )
        )
        # Partial NaN: should produce valid CI
        stats = reports.compute_gene_stats(bins, 0.15, interval_stats=["ci"])
        self.assertIn("ci_lo", stats)
        self.assertFalse(np.isnan(stats["ci_lo"]))

        # All NaN: should return empty stats (no crash)
        bins_allnan = bins.copy()
        bins_allnan["weight"] = np.nan
        stats2 = reports.compute_gene_stats(bins_allnan, 0.15, interval_stats=["ci"])
        self.assertNotIn("ci_lo", stats2)

    def test_squash_region_nan_gene(self):
        """squash_region handles NaN gene names without crashing (issue #900)."""
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1"],
                "start": [0, 1000],
                "end": [1000, 2000],
                "gene": ["GeneA", float("nan")],
                "log2": [0.1, 0.2],
                "probes": [5, 5],
                "weight": [1.0, 1.0],
            }
        )
        result = segfilters.squash_region(df)
        gene_val = result["gene"].iat[0]
        self.assertIsInstance(gene_val, str)
        self.assertNotIn("nan", gene_val.lower())
        self.assertEqual(gene_val, "GeneA")

        # All-NaN genes should produce the placeholder "-"
        df_allnan = df.copy()
        df_allnan["gene"] = [float("nan"), float("nan")]
        result2 = segfilters.squash_region(df_allnan)
        self.assertEqual(result2["gene"].iat[0], "-")

    def test_squash_by_groups_no_cross_region_merge(self):
        """squash_by_groups merges only contiguous runs (#677).

        The old code combined a run-id cumsum with an integer chromosome
        index by *addition*, which collided on position-unsorted input and
        merged non-consecutive bins -- yielding segments with start > end and
        a bogus end coordinate (the constant "10488" in the ticket).  This is
        the shared segment-construction path used by HMM, CBS, and all segment
        filters, so it is not specific to HMM or to the segmetrics command that
        first surfaced it.  Boundaries are now detected by adjacency.
        """
        # chr1 appears in two runs (300-400, then 100-200) with a chr2 bin
        # between them; all three share the same state/level.
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr2", "chr1"],
                "start": [300, 50, 100],
                "end": [400, 90, 200],
                "gene": ["-", "-", "-"],
                "log2": [0.0, 0.0, 0.0],
                "probes": [1, 1, 1],
                "weight": [1.0, 1.0, 1.0],
            }
        )
        cn = cnary.CopyNumArray(df)
        levels = pd.Series([2, 2, 2], index=cn.data.index)
        out = segfilters.squash_by_groups(cn, levels)
        # No segment may have start >= end ...
        self.assertTrue((out.start < out.end).all(), out.data)
        # ... and the two chr1 runs stay separate (3 segments, not 2) with all
        # probes conserved and no cross-chromosome absorption.
        self.assertEqual(len(out), 3)
        self.assertEqual(out["probes"].sum(), 3)

    def test_squash_region_unsorted_span(self):
        """squash_region spans [min start, max end] (#677).

        A run whose bins are not in ascending order must not yield start > end.
        """
        df = pd.DataFrame(
            {
                "chromosome": ["chr1", "chr1"],
                "start": [300, 100],
                "end": [400, 200],
                "gene": ["A", "A"],
                "log2": [0.0, 0.0],
                "probes": [1, 1],
                "weight": [1.0, 1.0],
            }
        )
        result = segfilters.squash_region(df)
        self.assertEqual(result["start"].iat[0], 100)
        self.assertEqual(result["end"].iat[0], 400)

    def test_squash_by_groups_cn1_cn2_nan(self):
        """Allele-specific squash splits on known cn1/cn2 changes but treats
        NaN (missing BAF) as 'no change', matching legacy grouping (#677).

        Without this, the adjacency check would split every NaN row off on its
        own (NaN != NaN), changing .cns output for `call`/`ampdel` on segments
        that lack BAF data (call sets cn1/cn2 = NaN there).
        """
        base = {
            "chromosome": ["chr1"] * 4,
            "start": [0, 100, 200, 300],
            "end": [100, 200, 300, 400],
            "gene": ["-"] * 4,
            "log2": [0.0] * 4,
            "probes": [1, 1, 1, 1],
            "weight": [1.0] * 4,
            "cn": [2, 2, 2, 2],
        }
        # A known cn1/cn2 change in the middle splits into 2 segments.
        known = cnary.CopyNumArray(
            pd.DataFrame({**base, "cn1": [1, 1, 2, 2], "cn2": [1, 1, 0, 0]})
        )
        out = segfilters.squash_by_groups(known, known["cn"])
        self.assertEqual(len(out), 2)
        # All-NaN cn1/cn2 (no BAF) stays merged into 1 segment, not 4.
        allnan = cnary.CopyNumArray(
            pd.DataFrame({**base, "cn1": [np.nan] * 4, "cn2": [np.nan] * 4})
        )
        out_nan = segfilters.squash_by_groups(allnan, allnan["cn"])
        self.assertEqual(len(out_nan), 1)
        self.assertTrue((out_nan.start < out_nan.end).all())

    def test_squash_by_groups_empty_and_single(self):
        """squash_by_groups handles empty and single-row input (#677)."""
        cols = ["chromosome", "start", "end", "gene", "log2", "probes", "weight"]
        empty = cnary.CopyNumArray(pd.DataFrame({c: [] for c in cols}))
        self.assertEqual(len(segfilters.squash_by_groups(empty, empty["log2"])), 0)
        one = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"],
                    "start": [10],
                    "end": [20],
                    "gene": ["G"],
                    "log2": [0.1],
                    "probes": [1],
                    "weight": [1.0],
                }
            )
        )
        out = segfilters.squash_by_groups(one, one["log2"])
        self.assertEqual(len(out), 1)
        self.assertEqual((out["start"].iat[0], out["end"].iat[0]), (10, 20))


class VATests(unittest.TestCase):
    """Tests for the VariantArray class."""

    def test_read(self):
        """Instantiate from a VCF file."""
        variants = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        self.assertGreater(len(variants), 0)

    def test_mirrored_baf_all_nan(self):
        """All-NaN frequencies mirror to all-NaN without warning (gh#407).

        np.nanmedian warns 'Mean of empty slice' on an all-NaN slice under
        older numpy; with pytest's filterwarnings=error that breaks the build
        (the 3.11-min job). _mirrored_baf must short-circuit before calling
        .median() on such a slice.
        """
        vals = pd.Series([np.nan, np.nan, np.nan])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = vary._mirrored_baf(vals)
        self.assertEqual(len(out), len(vals))
        self.assertTrue(out.isna().all())
        # A partially-NaN slice still mirrors normally (median over real values)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            mixed = vary._mirrored_baf(pd.Series([np.nan, 0.2, 0.3]))
        self.assertTrue(math.isnan(mixed.iloc[0]))
        # median 0.25 < 0.5 -> below-half branch (0.5 - shift); both already minor
        self.assertAlmostEqual(mixed.iloc[1], 0.2)
        self.assertAlmostEqual(mixed.iloc[2], 0.3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
