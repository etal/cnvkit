#!/usr/bin/env python
"""Tests for analysis commands: metrics, segmetrics, breaks, genemetrics, bintest, etc."""

import logging
import os
import shutil
import tempfile
import unittest
import warnings

import pytest

logging.basicConfig(level=logging.ERROR, format="%(message)s")
warnings.filterwarnings("ignore", category=DeprecationWarning)

import numpy as np
import pandas as pd
import pysam
from skgenome import tabio, GenomicArray as GA

import cnvlib
from conftest import linecount
from cnvlib import (
    access,
    antitarget,
    autobin,
    batch,
    bintest,
    call,
    cluster,
    cmdutil,
    cnary,
    commands,
    core,
    coverage,
    diagram,
    export,
    fix,
    heatmap,
    import_rna,
    importers,
    metrics,
    parallel,
    params,
    plots,
    purity,
    reference,
    reports,
    samutil,
    scatter,
    segfilters,
    segmentation,
    segmetrics,
    smoothing,
    vary,
)


class AnalysisTests(unittest.TestCase):
    """Tests for analysis commands: bintest, breaks, genemetrics, metrics, etc."""

    def test_bintest(self):
        """The 'bintest' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        # Simple
        rows = commands.do_bintest(cnarr, alpha=0.05)
        self.assertGreater(len(rows), 0)
        self.assertLess(len(rows), len(cnarr))
        # Versus segments
        rows = commands.do_bintest(cnarr, segarr, target_only=True)
        self.assertGreaterEqual(len(rows), len(segarr))
        self.assertLess(len(rows), len(cnarr))

    def test_breaks(self):
        """The 'breaks' command."""
        probes = cnvlib.read("formats/amplicon.cnr")
        segs = cnvlib.read("formats/amplicon.cns")
        rows = commands.do_breaks(probes, segs, 4)
        self.assertGreater(len(rows), 0)

    def test_genemetrics(self):
        """The 'genemetrics' command."""
        probes = cnvlib.read("formats/amplicon.cnr")
        rows = commands.do_genemetrics(probes, is_haploid_x_reference=True)
        self.assertGreater(len(rows), 0)
        segs = cnvlib.read("formats/amplicon.cns")
        rows = commands.do_genemetrics(
            probes, segs, 0.3, 4, is_haploid_x_reference=True
        )
        self.assertGreater(len(rows), 0)

    @pytest.mark.slow
    def test_genemetrics_with_stats(self):
        """The 'genemetrics' command with statistics options."""
        probes = cnvlib.read("formats/amplicon.cnr")
        segs = cnvlib.read("formats/amplicon.cns")

        # Test location statistics
        result = commands.do_genemetrics(
            probes,
            segs,
            0.3,
            4,
            is_haploid_x_reference=True,
            location_stats=["mean", "median", "mode", "p_ttest"],
        )
        self.assertGreater(len(result), 0)
        self.assertIn("mean", result.columns)
        self.assertIn("median", result.columns)
        self.assertIn("mode", result.columns)
        self.assertIn("p_ttest", result.columns)
        # P-values should be between 0 and 1 (excluding NaN)
        valid_p = result["p_ttest"].dropna()
        self.assertGreater(len(valid_p), 0)
        self.assertTrue((valid_p >= 0).all())
        self.assertTrue((valid_p <= 1).all())

        # Test spread statistics
        result = commands.do_genemetrics(
            probes,
            segs,
            0.3,
            4,
            is_haploid_x_reference=True,
            spread_stats=["stdev", "sem", "mad", "mse", "iqr", "bivar"],
        )
        self.assertGreater(len(result), 0)
        for stat in ["stdev", "sem", "mad", "mse", "iqr", "bivar"]:
            self.assertIn(stat, result.columns)
            # Spread statistics should be non-negative (excluding NaN)
            valid_stat = result[stat].dropna()
            if len(valid_stat) > 0:
                self.assertTrue((valid_stat >= 0).all())

        # Test interval statistics with adaptive smoothed bootstrap
        result = commands.do_genemetrics(
            probes,
            segs,
            0.3,
            4,
            is_haploid_x_reference=True,
            location_stats=["mean", "median"],
            interval_stats=["ci", "pi"],
            alpha=0.05,
            bootstraps=50,
            smoothed=10,  # Default threshold
        )
        self.assertGreater(len(result), 0)
        self.assertIn("ci_lo", result.columns)
        self.assertIn("ci_hi", result.columns)
        self.assertIn("pi_lo", result.columns)
        self.assertIn("pi_hi", result.columns)

        # Confidence intervals should be ordered: ci_lo <= mean <= ci_hi
        # Filter out rows with NaN
        valid_rows = result.dropna(subset=["ci_lo", "ci_hi", "mean"])
        self.assertGreater(len(valid_rows), 0)
        ci_contains_mean = (valid_rows["ci_lo"] <= valid_rows["mean"]) & (
            valid_rows["mean"] <= valid_rows["ci_hi"]
        )
        # Most should contain the mean (allow some misses due to random sampling)
        self.assertGreater(ci_contains_mean.sum(), len(valid_rows) * 0.9)

        # Prediction intervals should contain the mean
        valid_rows = result.dropna(subset=["pi_lo", "pi_hi", "mean"])
        self.assertGreater(len(valid_rows), 0)
        pi_contains_mean = (valid_rows["pi_lo"] <= valid_rows["mean"]) & (
            valid_rows["mean"] <= valid_rows["pi_hi"]
        )
        self.assertGreater(pi_contains_mean.sum(), len(valid_rows) * 0.9)

        # Intervals should be properly ordered (lo <= hi)
        valid_rows = result.dropna(subset=["pi_lo", "pi_hi", "ci_lo", "ci_hi"])
        self.assertGreater(len(valid_rows), 0)
        self.assertTrue((valid_rows["ci_lo"] <= valid_rows["ci_hi"]).all())
        self.assertTrue((valid_rows["pi_lo"] <= valid_rows["pi_hi"]).all())
        # PI is usually wider than CI, but may not always be for very small
        # segments with smoothed bootstrap, so just check most cases
        pi_width = valid_rows["pi_hi"] - valid_rows["pi_lo"]
        ci_width = valid_rows["ci_hi"] - valid_rows["ci_lo"]
        self.assertGreater((pi_width >= ci_width).sum(), len(valid_rows) * 0.8)

        # Test without segments (gene-level analysis)
        result = commands.do_genemetrics(
            probes,
            None,
            0.3,
            4,
            is_haploid_x_reference=True,
            location_stats=["mean", "median"],
            spread_stats=["stdev"],
            interval_stats=["pi"],
        )
        self.assertGreater(len(result), 0)
        self.assertIn("mean", result.columns)
        self.assertIn("stdev", result.columns)
        self.assertIn("pi_lo", result.columns)
        self.assertIn("pi_hi", result.columns)

        # Test different smoothed bootstrap thresholds
        # Small threshold (0) - should use BCa for all
        result_bca = commands.do_genemetrics(
            probes,
            segs,
            0.3,
            4,
            is_haploid_x_reference=True,
            location_stats=["mean"],
            interval_stats=["ci"],
            bootstraps=50,
            smoothed=0,  # Never smooth, always BCa
        )
        self.assertGreater(len(result_bca), 0)

        # Large threshold (1000) - should use smoothed bootstrap for all
        result_smooth = commands.do_genemetrics(
            probes,
            segs,
            0.3,
            4,
            is_haploid_x_reference=True,
            location_stats=["mean"],
            interval_stats=["ci"],
            bootstraps=50,
            smoothed=1000,  # Always smooth
        )
        self.assertGreater(len(result_smooth), 0)

        # CIs should be valid (lo <= hi)
        valid_bca = result_bca.dropna(subset=["ci_lo", "ci_hi"])
        valid_smooth = result_smooth.dropna(subset=["ci_lo", "ci_hi"])
        self.assertGreater(len(valid_bca), 0)
        self.assertGreater(len(valid_smooth), 0)
        self.assertTrue((valid_bca["ci_lo"] <= valid_bca["ci_hi"]).all())
        self.assertTrue((valid_smooth["ci_lo"] <= valid_smooth["ci_hi"]).all())

    def test_import_theta(self):
        """The 'import-theta' command."""
        cns = cnvlib.read("formats/nv3.cns")
        theta_fname = "formats/nv3.n3.results"
        for new_cns in commands.do_import_theta(cns, theta_fname):
            self.assertTrue(0 < len(new_cns) <= len(cns))

    def test_metrics(self):
        """The 'metrics' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segments = cnvlib.read("formats/amplicon.cns")
        result = metrics.do_metrics(cnarr, segments, skip_low=True)
        self.assertEqual(result.shape, (1, 6))
        values = result.loc[0, result.columns[1:]]
        for val in values:
            self.assertGreater(val, 0)

    def test_metrics_multisample_and_nan(self):
        """metrics over multiple samples, and with NaN bins, stays finite."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segments = cnvlib.read("formats/amplicon.cns")
        # Multiple samples -> one row each
        multi = metrics.do_metrics([cnarr, cnarr], [segments, segments], skip_low=True)
        self.assertEqual(multi.shape, (2, 6))
        # A NaN log2 bin must not poison the residual-scale estimators
        nan_cnarr = cnarr.copy()
        nan_cnarr.data.loc[nan_cnarr.data.index[0], "log2"] = np.nan
        res = metrics.do_metrics(nan_cnarr, segments, skip_low=True)
        scale_vals = res.loc[0, res.columns[2:]].to_numpy(dtype=float)
        self.assertTrue(np.all(np.isfinite(scale_vals)))

    def test_segmetrics_single_bin(self):
        """Single-probe segments get finite CIs bracketing the mean (no crash)."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        sm = segmetrics.do_segmetrics(
            cnarr,
            segarr,
            location_stats=["mean", "median"],
            spread_stats=["stdev"],
            interval_stats=["pi", "ci"],
            bootstraps=50,
            smoothed=True,
        )
        one = sm.data[sm.data["probes"] == 1]
        self.assertGreater(len(one), 0, "amplicon.cns has single-probe segments")
        self.assertTrue(
            ((one["ci_lo"] <= one["mean"]) & (one["mean"] <= one["ci_hi"])).all()
        )
        finite = np.isfinite(
            one[["mean", "ci_lo", "ci_hi", "pi_lo", "pi_hi"]].to_numpy(dtype=float)
        )
        self.assertTrue(finite.all())

    def test_segmetrics(self):
        """The 'segmetrics' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        sm = segmetrics.do_segmetrics(
            cnarr,
            segarr,
            location_stats=["mean", "median", "mode", "p_ttest"],
            spread_stats=["stdev", "sem", "iqr"],
            interval_stats=["pi", "ci"],
            bootstraps=50,
            smoothed=True,
        )
        # Restrict to segments with enough supporting probes for sane stats
        sm = sm[sm["probes"] > 3]
        self.assertTrue((sm["pi_lo"] < sm["median"]).all())
        self.assertTrue((sm["pi_hi"] > sm["median"]).all())
        self.assertTrue((sm["ci_lo"] < sm["mean"]).all())
        self.assertTrue((sm["ci_hi"] > sm["mean"]).all())

    @pytest.mark.slow
    def test_purity(self):
        """The 'purity' command."""
        segments = cnvlib.read("formats/tr95t.cns")

        # Basic call without VCF
        result = commands.do_purity(segments)
        self.assertIn("purity", result.columns)
        self.assertIn("ploidy", result.columns)
        self.assertIn("score", result.columns)
        self.assertGreater(len(result), 0)
        self.assertGreater(result["purity"].iloc[0], 0)
        self.assertLessEqual(result["purity"].iloc[0], 1.0)
        self.assertGreaterEqual(result["ploidy"].iloc[0], 1.5)
        self.assertLessEqual(result["ploidy"].iloc[0], 5.0)
        self.assertTrue((result["score"] > 0).all())

        # With VCF
        varr = tabio.read(
            "formats/na12878_na12882_mix.vcf",
            "vcf",
            skip_somatic=True,
        ).heterozygous()
        result_baf = commands.do_purity(segments, varr)
        self.assertEqual(list(result_baf.columns), ["purity", "ploidy", "score"])
        self.assertGreater(len(result_baf), 0)

        # Custom grid (coarser steps for speed)
        result_coarse = commands.do_purity(
            segments,
            min_purity=0.2,
            max_purity=0.8,
            purity_step=0.1,
            min_ploidy=2.0,
            max_ploidy=4.0,
            ploidy_step=0.5,
        )
        self.assertGreater(len(result_coarse), 0)
        self.assertGreaterEqual(result_coarse["purity"].iloc[0], 0.2)
        self.assertLessEqual(result_coarse["purity"].iloc[0], 0.8)

    def test_purity_file(self):
        """The 'purity' command: read_purity_tsv round-trip."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("purity\tploidy\tscore\n")
            f.write("0.65\t2.1\t42.5\n")
            f.write("0.70\t2.0\t41.0\n")
            tmp_path = f.name
        try:
            pur, plo = purity.read_purity_tsv(tmp_path)
            self.assertAlmostEqual(pur, 0.65)
            self.assertAlmostEqual(plo, 2.1)
        finally:
            os.unlink(tmp_path)
