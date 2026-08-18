#!/usr/bin/env python
"""Tests for analysis commands: metrics, segmetrics, breaks, genemetrics, bintest, etc."""

import inspect
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
from conftest import linecount

import cnvlib
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
from skgenome import GenomicArray as GA
from skgenome import tabio


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

    def test_breaks_names_one_gene_per_row(self):
        """A bin naming several genes yields one row per gene, not a pseudo-gene.

        The gene column is what the documented ``breaks ... | cut -f1 | sort -u``
        pipeline reads, so every value in it must be a name that exists as a
        gene. Tallying the whole comma-joined field invented five names on this
        fixture, among them 'GOPC,ROS1' -- under which a breakpoint in ROS1 was
        invisible to anyone looking for ROS1.
        """
        probes = cnvlib.read("formats/amplicon.cnr")
        segs = cnvlib.read("formats/amplicon.cns")
        rows = commands.do_breaks(probes, segs, 1)
        for gene in rows["gene"]:
            self.assertNotIn(",", gene)
        # The co-binned pair is reported once each, same breakpoint
        gopc = rows[rows["gene"] == "GOPC"]
        ros1 = rows[rows["gene"] == "ROS1"]
        self.assertEqual(len(gopc), 1)
        self.assertEqual(len(ros1), 1)
        for column in ("chromosome", "location", "change"):
            self.assertEqual(gopc[column].iat[0], ros1[column].iat[0])
        # ... and their probe counts differ, being their own loci's bins
        self.assertNotEqual(gopc["probes_left"].iat[0], ros1["probes_left"].iat[0])

    def test_breaks_probe_counts_match_genemetrics(self):
        """probes_left + probes_right is the locus's bin count in genemetrics.

        A copy number alteration affects a genomic region regardless of the gene
        model, so an off-target or placeholder bin inside a gene supports a
        breakpoint there as well as a targeted bin does. Both commands therefore
        count the whole run of bins assigned to the locus, and the sum is
        genemetrics' 'probes' exactly.
        """
        probes = cnvlib.read("formats/amplicon.cnr")
        segs = cnvlib.read("formats/amplicon.cns")
        rows = commands.do_breaks(probes, segs, 1)
        self.assertGreater(len(rows), 0)
        gene_probes: dict = {}
        genemetrics = commands.do_genemetrics(probes, threshold=-1e9, min_probes=0)
        for row in genemetrics.itertuples():
            gene_probes.setdefault((row.chromosome, row.gene), []).append(row.probes)
        for row in rows.itertuples():
            self.assertIn(
                row.probes_left + row.probes_right,
                gene_probes[row.chromosome, row.gene],
                f"{row.gene} at {row.chromosome}:{row.location}",
            )

    def test_breaks_ignores_gap_between_two_loci(self):
        """A segment boundary between two loci of one name is not a breakpoint.

        Condensing every occurrence of a name on a chromosome into one interval
        put the gap between two loci inside it, so a boundary there was reported
        as a breakpoint 'within the gene', with the bins of one locus counted on
        the left and the other's on the right.
        """
        probes = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 1000, 2000, "REPEAT", 0.0],
                ["chr1", 2000, 3000, "REPEAT", 0.0],
                ["chr1", 3000, 4000, "MIDGENE", 0.0],
                ["chr1", 4000, 5000, "MIDGENE", 0.0],
                ["chr1", 9000, 9500, "REPEAT", 0.0],
                ["chr1", 9500, 10000, "REPEAT", 0.0],
            ]
        )
        segs = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 1000, 3500, "-", 0.0],
                ["chr1", 3500, 10000, "-", 1.0],
            ]
        )
        rows = commands.do_breaks(probes, segs, 1)
        # The boundary at 3500 lies inside MIDGENE, and inside nothing else
        self.assertEqual(list(rows["gene"]), ["MIDGENE"])

    def test_breaks_counts_placeholder_bins_within_a_gene(self):
        """Placeholder bins inside a gene count toward its breakpoint support.

        A bin labeled 'CGH' is a backbone bait, not a gene, and one landing
        inside a real gene is still a measurement of that locus. So it neither
        splits the gene nor goes uncounted: here two of the three bins left of
        the breakpoint carry the gene's name and the third is the placeholder.
        """
        probes = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 1000, 2000, "GENE", 0.0],
                ["chr1", 2000, 3000, "CGH", 0.0],
                ["chr1", 3000, 4000, "GENE", 0.0],
                ["chr1", 4000, 5000, "GENE", 0.0],
            ]
        )
        segs = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 1000, 4000, "-", 0.0],
                ["chr1", 4000, 5000, "-", 1.0],
            ]
        )
        rows = commands.do_breaks(probes, segs, 1)
        self.assertEqual(list(rows["gene"]), ["GENE"])
        self.assertEqual(rows["probes_left"].iat[0], 3)
        self.assertEqual(rows["probes_right"].iat[0], 1)

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
        """Single-probe segments get finite CIs bracketing the mean (no crash).

        Passes ``smoothed=True`` deliberately: no command line can produce a bool
        any more, so the bool arm of the ``bool | int`` union is reachable only
        from the test suite. The single-probe rows asserted on here never reach
        the smoother -- ``confidence_interval_bootstrap`` returns a zero-width
        interval before any resampling -- so what this pins is that the early
        return survives the smoothed configuration.
        """
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

    def test_smoothed_bootstrap_ci_is_reproducible(self):
        """Repeated ``--ci`` runs on the same input must give identical bounds.

        Both draws inside ``confidence_interval_bootstrap`` -- the resample
        indices and the Gaussian smoothing noise -- are taken from one seeded
        generator, so repeated runs on identical input agree exactly. The noise
        draw was unseeded between the NumPy modernization and this test, which
        gave ``segmetrics --ci`` and ``genemetrics --ci`` run-to-run variation in
        ``ci_lo``/``ci_hi`` for every segment at or below the
        ``--smooth-bootstrap`` threshold -- 41 of amplicon's 80 segments, enough
        to move a ``call --filter ci`` decision on a handful of them per run.

        The weights are deliberately below 1. The noise scales as ``sqrt(1 - w)``,
        so unit weights annihilate it and the unseeded code passed a
        reproducibility assertion written that way.
        """
        values = np.array([-0.3, 0.1, 0.25])
        weights = np.array([0.6, 0.8, 0.45])
        ci_args = (values, weights, 0.05)
        first = segmetrics.confidence_interval_bootstrap(*ci_args, smoothed=True)
        self.assertGreater(first[1] - first[0], 0, "smoothing must widen the interval")
        second = segmetrics.confidence_interval_bootstrap(*ci_args, smoothed=True)
        np.testing.assert_array_equal(first, second)
        # Same guarantee through the command, where the threshold rather than the
        # bool selects the smoothed path.
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        n_smoothed = ((segarr["probes"] >= 2) & (segarr["probes"] <= 10)).sum()
        self.assertGreater(n_smoothed, 40, "fixture must exercise the smoothed path")
        runs = [
            segmetrics.do_segmetrics(cnarr, segarr, interval_stats=["ci"])
            for _ in range(2)
        ]
        for col in ("ci_lo", "ci_hi"):
            np.testing.assert_array_equal(runs[0][col], runs[1][col])

    def test_same_size_segments_do_not_share_a_resample_pattern(self):
        """Two segments of equal bin count must not draw the same indices.

        One fixed seed per call made the resample index matrix a function of
        ``(n_bins, bootstraps)`` alone, which left the unsmoothed interval
        equivariant under an increasing affine map of the values: the weighted
        resample means commute with ``a*x + b``, and so do BCa's jackknife
        influence values, its bias correction (a count of means below the
        observed one) and its scale-free acceleration. So
        ``ci(2v + 0.5) == 2*ci(v) + 0.5`` held to machine precision -- possible
        only while both calls share a stream. Measured on the pre-fix code, the
        two sides agreed to better than 1e-12.

        This is the property to defend rather than "the bounds differ", because
        it fails for any rewrite that keys the stream on the bin count, the
        bootstrap count, or anything else two same-size segments have in common.
        """
        rng = np.random.default_rng(7)
        values = rng.standard_normal(20) * 0.3
        weights = rng.uniform(0.4, 0.9, 20)
        scale, shift = 2.0, 0.5
        plain = segmetrics.confidence_interval_bootstrap(values, weights, 0.05)
        mapped = segmetrics.confidence_interval_bootstrap(
            scale * values + shift, weights, 0.05
        )
        self.assertFalse(
            np.allclose(mapped, scale * plain + shift, rtol=1e-6, atol=1e-6),
            f"{mapped} is the affine image of {plain}, so both segments drew the "
            "same bootstrap indices",
        )
        # Both intervals are still intervals around their own weighted mean.
        for ci, vals in ((plain, values), (mapped, scale * values + shift)):
            mean = np.average(vals, weights=weights)
            self.assertLess(ci[0], mean)
            self.assertLess(mean, ci[1])

    def test_ci_does_not_depend_on_the_rest_of_the_file(self):
        """``--ci`` on one chromosome reproduces that chromosome's rows exactly.

        The seed is derived from each segment's own values, so a segment's bounds
        cannot depend on its position in the input or on what else the input
        contained. Decorrelating same-size segments by counting iterations would
        have passed the reproducibility test above and broken this one, which is
        why it is asserted rather than left as a measured property.
        """
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        whole = segmetrics.do_segmetrics(cnarr, segarr, interval_stats=["ci"])
        part = segmetrics.do_segmetrics(
            cnarr.in_range("chr7"), segarr.in_range("chr7"), interval_stats=["ci"]
        )
        self.assertGreater(len(part), 1, "fixture must have several chr7 segments")
        whole_chr7 = whole.in_range("chr7")
        for col in ("ci_lo", "ci_hi"):
            np.testing.assert_array_equal(whole_chr7[col], part[col])

    def test_ci_ignores_the_weight_columns_dtype(self):
        """Weights read as integers must give the bounds their float values give.

        Reachable through CNVkit's own writer, not only a hand-written file: it
        formats with ``%.6g``, so an all-1.0 weight column is written as ``1`` and
        reads back as int64. ``_rng_from_values`` casts before digesting for this
        reason; the consequences of not casting are recorded there.
        """
        values = np.array([-0.3, 0.1, 0.25, 0.4, -0.15])
        as_int = np.ones(len(values), dtype=np.int64)
        as_float = as_int.astype(np.float64)
        self.assertNotEqual(as_int.dtype, as_float.dtype, "case must differ in dtype")
        np.testing.assert_array_equal(
            segmetrics.confidence_interval_bootstrap(values, as_int, 0.05),
            segmetrics.confidence_interval_bootstrap(values, as_float, 0.05),
        )

    def test_bca_jackknife_matches_leave_one_out(self):
        """BCa's influence values are the leave-one-out weighted means.

        The closed form over the totals agrees with the definition only if the
        weights are carried through, so this pins the algebra rather than the
        speed: a version that forgot the weights would still be fast.
        """
        cnarr = cnvlib.read("formats/amplicon.cnr")
        values = cnarr["log2"].to_numpy()[:200]
        weights = cnarr["weight"].to_numpy()[:200]
        self.assertGreater(weights.std(), 0, "constant weights would not discriminate")
        expected = np.array(
            [
                np.average(np.delete(values, i), weights=np.delete(weights, i))
                for i in range(len(values))
            ]
        )
        observed = segmetrics._jackknife_weighted_means(values, weights)
        np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=0)

    def test_smoothing_defaults_agree_across_the_relay(self):
        """Every layer that carries ``smoothed`` must default to the same value.

        The commit that turned ``--smooth-bootstrap`` from a flag into a bin-count
        threshold moved six signature defaults and missed
        ``confidence_interval_bootstrap``'s, which stayed at ``False``. That miss
        is invisible at the call sites, since each layer passes the argument on
        explicitly; it would surface only for a direct API user, as segments that
        silently stopped being smoothed. ``bootstraps`` rides the same relay and
        is checked here for the same reason.
        """
        relay = (
            segmetrics.do_segmetrics,
            segmetrics.confidence_interval_bootstrap,
            reports.do_genemetrics,
            reports.gene_metrics_by_gene,
            reports.gene_metrics_by_segment,
            reports.compute_gene_stats,
            reports.group_by_genes,
        )
        for param, dest, expected in (
            ("smoothed", "smooth_bootstrap", 10),
            ("bootstraps", "bootstrap", 100),
        ):
            for func in relay:
                default = inspect.signature(func).parameters[param].default
                self.assertEqual(
                    default,
                    expected,
                    f"{func.__module__}.{func.__name__} defaults {param} to "
                    f"{default!r}, not {expected!r}",
                )
            cli_defaults = {
                action.default
                for parser in (commands.P_segmetrics, commands.P_genemetrics)
                for action in parser._actions
                if action.dest == dest
            }
            self.assertEqual(
                cli_defaults,
                {expected},
                f"--{dest.replace('_', '-')} must default to {expected} at both "
                "subcommands that offer it, matching the library signatures above "
                "and the numbers quoted in doc/reports.rst",
            )

    def test_segmetrics(self):
        """The 'segmetrics' command, at the smoothing default users get."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        sm = segmetrics.do_segmetrics(
            cnarr,
            segarr,
            location_stats=["mean", "median", "mode", "p_ttest"],
            spread_stats=["stdev", "sem", "iqr"],
            interval_stats=["pi", "ci"],
            bootstraps=50,
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
