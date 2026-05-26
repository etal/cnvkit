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
from cnvlib.cnary import CopyNumArray
from cnvlib import (
    cmdutil,
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

    def test_guess_xx_indeterminate_defaults_female(self):
        """Indeterminate sex (no chrX) defaults to female, never silently male.

        `guess_xx` still returns None -- an honest "can't tell" that the
        reference target/antitarget reconciliation relies on -- but the
        decision-level helper `is_female_default` and `verify_sample_sex`
        collapse None to female. Female is the safe default: guessing male
        for a sample whose X is actually diploid would spuriously inflate
        chrX by +1 (the gh#360 failure family).
        """
        from cnvlib.cnary import is_female_default

        # The helper: None -> female; real guesses pass through as plain bool.
        self.assertIs(is_female_default(None), True)
        self.assertIs(is_female_default(np.bool_(True)), True)
        self.assertIs(is_female_default(np.bool_(False)), False)

        # An autosomes-only sample: guess_xx genuinely can't tell (None) ...
        auto_only = CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 5 + ["chr2"] * 5,
                    "start": list(range(0, 5000, 1000)) * 2,
                    "end": list(range(1000, 6000, 1000)) * 2,
                    "gene": ["-"] * 10,
                    "log2": [0.0] * 10,
                    "weight": [1.0] * 10,
                }
            )
        )
        self.assertIsNone(auto_only.guess_xx(verbose=False))
        # ... but verify_sample_sex returns female (a real bool), not None/male.
        result = cmdutil.verify_sample_sex(auto_only, None, False, None)
        self.assertIs(result, True)
        # An explicit override still wins, with no misleading mismatch warning.
        self.assertIs(cmdutil.verify_sample_sex(auto_only, "male", False, None), False)

    def test_guess_xx_y_presence_required_for_male(self):
        """Without positive chrY evidence, even a haploid-looking chrX defers to female (#954).

        WGS with ``--diploid-parx-genome grch38`` strips chrY to NULL_LOG2_COVERAGE
        for every sample, so chrY carries no signal. The published log in #954
        shows female samples whose chrX ratios drifted into the male hemisphere
        (~ -0.5 to -1.2) getting mis-called male by the old multiplicative
        combination of chrX and chrY maleness ratios. Per the design principle
        "Y-presence is the best proxy for the germline baseline", a sample
        with no detectable chrY should not be called male on chrX noise alone.

        The boundary cases below mirror the actual #954 ratios; ground truth for
        all of them is female.
        """

        def _make(chrx_ratio, chry_ratio):
            # 200 autosome bins around log2=0 (anchors the autosome median),
            # 50 chrX bins around chrx_ratio, 50 chrY bins around chry_ratio.
            rng = np.random.default_rng(0)
            rows = []
            for c in ("chr1", "chr2"):
                for i in range(200):
                    rows.append(
                        (
                            c,
                            i * 1000,
                            i * 1000 + 100,
                            "-",
                            float(rng.normal(0, 0.1)),
                            1.0,
                        )
                    )
            for i in range(50):
                rows.append(
                    (
                        "chrX",
                        i * 1000,
                        i * 1000 + 100,
                        "-",
                        float(rng.normal(chrx_ratio, 0.1)),
                        1.0,
                    )
                )
            for i in range(50):
                rows.append(
                    (
                        "chrY",
                        i * 1000,
                        i * 1000 + 100,
                        "-",
                        float(rng.normal(chry_ratio, 0.1)),
                        1.0,
                    )
                )
            return CopyNumArray(
                pd.DataFrame(
                    rows,
                    columns=["chromosome", "start", "end", "gene", "log2", "weight"],
                )
            )

        # #954's mis-called samples: chrX between -0.5 and -1.2, chrY effectively
        # NULL (-21 to -23). All should now be female.
        for chrx, chry in [
            (-0.501, -22.0),  # the knife-edge sample in #954
            (-0.843, -23.1),
            (-0.838, -22.0),
            (-1.05, -21.3),
            (-1.21, -23.3),
        ]:
            with self.subTest(chrx=chrx, chry=chry):
                cna = _make(chrx, chry)
                self.assertIs(
                    cna.guess_xx(verbose=False),
                    np.True_,
                    f"chrX={chrx} chrY={chry} should be female (no Y evidence)",
                )

        # Sanity: a sample with a haploid-looking chrX AND real chrY signal
        # (~ -1, i.e. half-coverage, not the NULL floor) is still confidently male.
        cna = _make(chrx_ratio=-1.0, chry_ratio=-1.0)
        self.assertIs(cna.guess_xx(verbose=False), np.False_)

    def test_guess_xx_consistent_between_cnr_and_cns(self):
        """Same chrX/chrY medians must give the same sex call regardless of
        bin/segment count (#785).

        The pre-existing Mood's-median-test maleness ratios depended on the
        chi-square magnitude, which scales with cell counts, so a .cnr with
        thousands of bins and a .cns with a handful of segments could produce
        different chrx_male_lr / chry_male_lr for identical chrx_ratio /
        chry_ratio. The ratio-based maleness uses only the median differences,
        so identical summaries -> identical sex call.
        """
        rng = np.random.default_rng(0)

        def _build(n_auto, n_x, n_y, *, x_log2=-1.03, y_log2=0.236):
            rows = []
            for c in ("chr1", "chr2"):
                for i in range(n_auto):
                    rows.append(
                        (
                            c,
                            i * 1000,
                            i * 1000 + 100,
                            "-",
                            float(rng.normal(0, 0.05)),
                            1.0,
                        )
                    )
            for i in range(n_x):
                rows.append(
                    (
                        "chrX",
                        i * 1000,
                        i * 1000 + 100,
                        "-",
                        float(rng.normal(x_log2, 0.05)),
                        1.0,
                    )
                )
            for i in range(n_y):
                rows.append(
                    (
                        "chrY",
                        i * 1000,
                        i * 1000 + 100,
                        "-",
                        float(rng.normal(y_log2, 0.05)),
                        1.0,
                    )
                )
            return CopyNumArray(
                pd.DataFrame(
                    rows,
                    columns=["chromosome", "start", "end", "gene", "log2", "weight"],
                )
            )

        # The "many-bin" view (think .cnr) and the "few-segment" view (think
        # .cns) of the same sample. With X log2 ~ -1.03 and Y log2 ~ +0.236,
        # this is the clearly-male case from #785 -- both representations must
        # agree.
        cnr_like = _build(n_auto=500, n_x=200, n_y=80)
        cns_like = _build(n_auto=20, n_x=5, n_y=3)
        self.assertIs(cnr_like.guess_xx(verbose=False), np.False_)  # male
        self.assertIs(cns_like.guess_xx(verbose=False), np.False_)  # male
        self.assertEqual(
            cnr_like.guess_xx(verbose=False),
            cns_like.guess_xx(verbose=False),
            ".cnr-like (many bins) and .cns-like (few segments) with same "
            "chrX/chrY medians must agree on sex (#785).",
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


class CNAryByGeneTests(unittest.TestCase):
    """Tests for the by_gene method with various index types."""

    def test_by_gene_with_numeric_index(self):
        """Test by_gene with standard numeric index (default)."""
        # Create test data with numeric index
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 6,
                "start": [100, 200, 300, 400, 500, 600],
                "end": [150, 250, 350, 450, 550, 650],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2, -0.3],
            }
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene and check for duplicates
        gene_groups = list(cnarr.by_gene())
        all_coords = []

        for gene_name, gene_cnarr in gene_groups:
            for row in gene_cnarr.data.itertuples():
                coord = (row.chromosome, row.start, row.end)
                self.assertNotIn(
                    coord,
                    all_coords,
                    f"Duplicate coordinate found in {gene_name}: {coord}",
                )
                all_coords.append(coord)

        # Should have processed all 6 rows
        self.assertEqual(len(all_coords), 6)

    def test_by_gene_with_non_sequential_index(self):
        """Test by_gene with non-sequential numeric index."""
        # Create test data with gaps in the index
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [100, 200, 300, 400, 500],
                "end": [150, 250, 350, 450, 550],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2],
            },
            index=[0, 5, 10, 15, 20],  # Non-sequential index
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene and check for duplicates
        gene_groups = list(cnarr.by_gene())
        all_coords = []

        for gene_name, gene_cnarr in gene_groups:
            for row in gene_cnarr.data.itertuples():
                coord = (row.chromosome, row.start, row.end)
                self.assertNotIn(
                    coord,
                    all_coords,
                    f"Duplicate coordinate found in {gene_name}: {coord}",
                )
                all_coords.append(coord)

        # Should have processed all 5 rows
        self.assertEqual(len(all_coords), 5)

    def test_by_gene_preserves_boundaries(self):
        """Test that by_gene correctly handles gene boundaries without overlap."""
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [100, 200, 300, 400, 500],
                "end": [150, 250, 350, 450, 550],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2],
            }
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        gene_groups = list(cnarr.by_gene())

        # Check that we have the expected groups
        gene_names = [name for name, _ in gene_groups]
        self.assertIn("GeneA", gene_names)
        self.assertIn("GeneB", gene_names)
        self.assertIn("Antitarget", gene_names)

        # Check GeneA has exactly 2 rows
        gene_a = next(ga for name, ga in gene_groups if name == "GeneA")
        self.assertEqual(len(gene_a), 2)
        self.assertEqual(gene_a.data.iloc[0]["start"], 100)
        self.assertEqual(gene_a.data.iloc[1]["end"], 250)

        # Check GeneB has exactly 2 rows
        gene_b = next(ga for name, ga in gene_groups if name == "GeneB")
        self.assertEqual(len(gene_b), 2)
        self.assertEqual(gene_b.data.iloc[0]["start"], 400)
        self.assertEqual(gene_b.data.iloc[1]["end"], 550)

    def test_by_gene_with_duplicate_index_labels(self):
        """Test by_gene with duplicate index labels (edge case for get_loc)."""
        # Create test data with duplicate index labels
        # This simulates the case where get_loc() returns a slice or boolean array
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 6,
                "start": [100, 200, 300, 400, 500, 600],
                "end": [150, 250, 350, 450, 550, 650],
                "gene": ["GeneA", "GeneA", "Antitarget", "GeneB", "GeneB", "GeneB"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2, -0.3],
            },
            index=[0, 1, 1, 2, 3, 3],  # Duplicate labels at positions 1-2 and 5-6
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene and verify no duplicates
        gene_groups = list(cnarr.by_gene())
        all_coords = []

        for gene_name, gene_cnarr in gene_groups:
            for row in gene_cnarr.data.itertuples():
                coord = (row.chromosome, row.start, row.end)
                self.assertNotIn(
                    coord,
                    all_coords,
                    f"Duplicate coordinate found in {gene_name}: {coord}",
                )
                all_coords.append(coord)

        # Should have processed all 6 rows
        self.assertEqual(len(all_coords), 6)

    def test_by_gene_multiple_genes_per_bin(self):
        """Test by_gene with bins spanning multiple genes (comma-separated)."""
        data = pd.DataFrame(
            {
                "chromosome": ["chr1"] * 5,
                "start": [100, 200, 300, 400, 500],
                "end": [150, 250, 350, 450, 550],
                "gene": ["GeneA", "GeneA,GeneB", "GeneB", "Antitarget", "GeneC"],
                "log2": [0.1, 0.2, 0.0, -0.1, -0.2],
            }
        )
        cnarr = CopyNumArray(data, {"sample_id": "test"})

        # Group by gene
        gene_groups = list(cnarr.by_gene())
        gene_dict = {name: ga for name, ga in gene_groups}

        # The bin at position 1 with "GeneA,GeneB" should appear in both genes
        self.assertIn("GeneA", gene_dict)
        self.assertIn("GeneB", gene_dict)
        self.assertIn("GeneC", gene_dict)

        # GeneA should include the shared bin
        gene_a_coords = [
            (row.chromosome, row.start, row.end)
            for row in gene_dict["GeneA"].data.itertuples()
        ]
        self.assertIn(("chr1", 200, 250), gene_a_coords)

        # GeneB should also include the shared bin
        gene_b_coords = [
            (row.chromosome, row.start, row.end)
            for row in gene_dict["GeneB"].data.itertuples()
        ]
        self.assertIn(("chr1", 200, 250), gene_b_coords)


if __name__ == "__main__":
    unittest.main()


if __name__ == "__main__":
    unittest.main(verbosity=2)
