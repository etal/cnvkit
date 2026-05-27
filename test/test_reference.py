#!/usr/bin/env python
"""Tests for the reference command (and clustering used to build references)."""

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


def _write_coverage_cnn(
    path, chrx_log2, chry_log2, chrx_sd, chry_sd, n_auto=200, n_x=40, n_y=40, seed=7
):
    """Write a synthetic coverage .cnn with tunable sex-chromosome signal.

    Used to reproduce target/antitarget sex-inference conflicts: autosomes sit
    at log2 0, and chrX/chrY are drawn around the given means so a caller can
    make, e.g., a haploid chrX with a near-empty (deeply negative) chrY.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for chrom, n in [("chr1", n_auto), ("chr2", n_auto)]:
        for i in range(n):
            rows.append(
                (
                    chrom,
                    i * 1000,
                    i * 1000 + 1000,
                    "Background",
                    100,
                    float(rng.normal(0, 0.2)),
                )
            )
    for i in range(n_x):
        rows.append(
            (
                "chrX",
                i * 1000,
                i * 1000 + 1000,
                "Background",
                50,
                float(rng.normal(chrx_log2, chrx_sd)),
            )
        )
    for i in range(n_y):
        rows.append(
            (
                "chrY",
                i * 1000,
                i * 1000 + 1000,
                "Background",
                2,
                float(rng.normal(chry_log2, chry_sd)),
            )
        )
    pd.DataFrame(
        rows, columns=["chromosome", "start", "end", "gene", "depth", "log2"]
    ).to_csv(path, sep="\t", index=False)


class ReferenceTests(unittest.TestCase):
    """Tests for the reference and fix commands."""

    def test_reference(self):
        """The 'reference' command."""
        # Empty/unspecified antitargets
        nlines = linecount("formats/amplicon.cnr") - 1
        ref = commands.do_reference(["formats/amplicon.cnr"], ["formats/empty"])
        self.assertEqual(len(ref), nlines)
        ref = commands.do_reference(["formats/amplicon.cnr"])
        self.assertEqual(len(ref), nlines)
        # Empty/unspecified antitargets, flat reference
        nlines = linecount("formats/amplicon.bed")
        ref = commands.do_reference_flat("formats/amplicon.bed", "formats/empty")
        self.assertEqual(len(ref), nlines)
        ref = commands.do_reference_flat("formats/amplicon.bed")
        self.assertEqual(len(ref), nlines)
        # Misc
        ref = cnvlib.read("formats/reference-tr.cnn")
        targets, antitargets = reference.reference2regions(ref)
        self.assertLess(0, len(antitargets))
        self.assertEqual(len(antitargets), (ref["gene"] == "Background").sum())
        self.assertEqual(len(targets), len(ref) - len(antitargets))

    def test_reference_gender_input(self):
        """Test whether correct log2-ratios are calculated for sex chromosomes in reference command"""
        expected_haploid_log2 = [-1.0]
        expected_diploid_log2 = [0.0]

        # Test case when sample gender is provided by user
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"], female_samples=False
        )  # is_haploid_x_reference=False
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"], female_samples=True
        )  # is_haploid_x_reference=False
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_diploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_diploid_log2,
        )
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"],
            female_samples=False,
            is_haploid_x_reference=True,
        )
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"],
            female_samples=True,
            is_haploid_x_reference=True,
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_haploid_log2,
        )

        # Test case when sample gender is guessed from input
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"]
        )  # is_haploid_x_reference=False
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"]
        )  # is_haploid_x_reference=False
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_diploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_diploid_log2,
        )
        ref_male = commands.do_reference(
            ["formats/ref_test_male.cnn"], is_haploid_x_reference=True
        )
        ref_female = commands.do_reference(
            ["formats/ref_test_female.cnn"], is_haploid_x_reference=True
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_x_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_male["log2"][ref_male.chr_y_filter()].to_list(), expected_haploid_log2
        )
        self.assertListEqual(
            ref_female["log2"][ref_female.chr_x_filter()].to_list(),
            expected_haploid_log2,
        )

    def test_shift_sex_chroms_unknown_defaults_female(self):
        """A sample absent from the inferred-sexes dict is treated as female.

        `infer_sexes` deliberately omits samples whose sex couldn't be
        determined (preserving a None signal for target/antitarget
        reconciliation), so `shift_sex_chroms` must apply the safe default
        itself: no +1 shift on chrX. The old code fell through to the male
        branch and silently inflated chrX (Defect A / #360 family).
        """
        cnarr = cnvlib.read("formats/ref_test_female.cnn")
        is_chr_x = cnarr.chr_x_filter()
        is_chr_y = cnarr.chr_y_filter()
        before_x = cnarr["log2"][is_chr_x].to_numpy().copy()
        # Empty sexes dict => this sample's sex is unknown.
        reference.shift_sex_chroms(cnarr, {}, np.zeros(len(cnarr)), is_chr_x, is_chr_y)
        after_x = cnarr["log2"][is_chr_x].to_numpy()
        # Female default: chrX unchanged (would be +1 under the None->male bug).
        np.testing.assert_allclose(after_x, before_x)

    def test_reconcile_sex_guesses(self):
        """Target/antitarget sex conflicts resolve by chrX confidence (#846).

        chrX is well-sampled in both targets and antitargets, whereas
        antitarget chrY is frequently too sparse to carry signal -- an empty
        antitarget chrY can drag a real male call to female. So when the two
        sources disagree, trust whichever source's chrX 'maleness' ratio is
        furthest from the male/female boundary (1.0). Each guess is
        ``(is_xx, chrx_male_lr)``; ``chrx_male_lr > 1`` favors male.
        """
        reconcile = reference._reconcile_sex_guesses
        MALE, FEMALE = np.False_, np.True_  # is_xx: True == female

        # The #846 bug: target clearly male, antitarget flipped female by an
        # empty chrY (its chrX is far less decisive). Keep the target's male call.
        out = reconcile({"S1": (MALE, 499.0)}, {"S1": (FEMALE, 0.29)})
        self.assertFalse(out["S1"])

        # The legitimate rescue this reconciliation must preserve: target chrX
        # was too sparse to look male, antitarget chrX decisively haploid.
        out = reconcile({"S2": (FEMALE, 0.9)}, {"S2": (MALE, 5.0)})
        self.assertFalse(out["S2"])

        # Agreement is untouched; confidence is irrelevant.
        out = reconcile({"S3": (FEMALE, 0.5)}, {"S3": (FEMALE, 0.6)})
        self.assertTrue(out["S3"])

        # Target couldn't tell (absent) -> fall back to the antitarget guess.
        out = reconcile({}, {"S4": (MALE, 3.0)})
        self.assertFalse(out["S4"])

        # A tie favors the target (its sex chromosomes are deliberately baited).
        out = reconcile({"S5": (FEMALE, 0.5)}, {"S5": (MALE, 2.0)})
        self.assertTrue(out["S5"])

        # When the target wins, its FULL chrX+chrY call stands -- we must not
        # re-derive sex from chrX alone. Here the target's chrX is weak (0.8,
        # chrX-only would say female) but its reliable chrY made the combined
        # call male; the sparse antitarget must not strip that chrY evidence.
        out = reconcile({"S6": (MALE, 0.8)}, {"S6": (FEMALE, 0.9)})
        self.assertFalse(out["S6"])

        # When the antitarget wins, only its chrX counts (its chrY is junk): both
        # sources' chrX read male, but the antitarget's sparse chrY flipped its
        # combined call to female. Trusting that combined call would re-break #846.
        out = reconcile({"S7": (MALE, 1.1)}, {"S7": (FEMALE, 1.5)})
        self.assertFalse(out["S7"])

    def test_reference_antitarget_empty_chry_keeps_male(self):
        """do_reference's inference must not flip a male sample to female when
        the antitarget chrY is near-empty (#846, #863).

        Builds a clearly-male target (haploid chrX, covered chrY) and an
        antitarget whose chrY is near-absent -- the configuration that made
        antitarget inference report 'female' and override the target. The fix
        keeps the male call without users having to delete antitarget chrY.
        """
        tmpdir = tempfile.mkdtemp()
        try:
            tgt = os.path.join(tmpdir, "S1.targetcoverage.cnn")
            anti = os.path.join(tmpdir, "S1.antitargetcoverage.cnn")
            _write_coverage_cnn(
                tgt,
                chrx_log2=-1.0,
                chry_log2=-1.0,
                chrx_sd=0.2,
                chry_sd=0.3,
                n_x=60,
                n_y=40,
            )
            # Antitarget chrY: deep below the NULL_LOG2_COVERAGE-leaning floor
            # for "absent" (well past the female/male midpoint at log2=-10), so
            # the chrY maleness gate flips the antitarget call to female even
            # though the antitarget chrX still reads haploid. This is the
            # configuration that fooled the pre-#1087 reconciliation into
            # preferring the antitarget (#846).
            _write_coverage_cnn(
                anti,
                chrx_log2=-0.5,
                chry_log2=-15.0,
                chrx_sd=0.6,
                chry_sd=1.5,
                n_x=40,
                n_y=80,
            )

            t_guesses = reference._guess_sexes([tgt], False, None, verbose=False)
            a_guesses = reference._guess_sexes([anti], False, None, verbose=False)
            # Precondition: this fixture really does trigger the conflict.
            self.assertFalse(t_guesses["S1"][0])  # target -> male
            self.assertTrue(a_guesses["S1"][0])  # antitarget -> female (the trap)

            sexes = reference._reconcile_sex_guesses(t_guesses, a_guesses)
            self.assertFalse(sexes["S1"])  # reconciled -> male, not flipped
        finally:
            shutil.rmtree(tmpdir)

    def test_fix(self):
        """The 'fix' command."""
        # Extract fake target/antitarget bins from a combined file
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        rng = np.random.default_rng(42)  # Use fixed seed for reproducible tests
        tgt_bins.log2 += rng.standard_normal(len(tgt_bins)) / 5
        anti_bins = ref[is_bg]
        anti_bins.log2 += rng.standard_normal(len(anti_bins)) / 5
        blank_bins = cnary.CopyNumArray([])
        # Typical usage (hybrid capture)
        cnr = commands.do_fix(tgt_bins, anti_bins, ref)
        self.assertTrue(0 < len(cnr) <= len(ref))
        # Blank antitargets (WGS or amplicon)
        cnr = commands.do_fix(tgt_bins, blank_bins, ref[~is_bg])
        self.assertTrue(0 < len(cnr) <= len(tgt_bins))

    def test_fix_bias_smoother_loess(self):
        """do_fix accepts bias_smoother='loess' end-to-end (gh#1028).

        Exercises the full plumb-through: do_fix -> load_adjust_coverages ->
        center_by_window -> smoothing.loess. The LOESS path must produce a
        finite .cnr of the same length as the default median path, and the
        log2 values must actually differ -- otherwise the parameter is being
        silently dropped somewhere in the stack.
        """
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        rng = np.random.default_rng(1028)
        tgt_bins.log2 += rng.standard_normal(len(tgt_bins)) / 5
        anti_bins = ref[is_bg]
        anti_bins.log2 += rng.standard_normal(len(anti_bins)) / 5

        cnr_median = commands.do_fix(tgt_bins, anti_bins, ref, bias_smoother="median")
        cnr_loess = commands.do_fix(tgt_bins, anti_bins, ref, bias_smoother="loess")

        # Both paths produce same-length, finite output
        self.assertEqual(len(cnr_loess), len(cnr_median))
        self.assertTrue(np.isfinite(cnr_loess["log2"]).all())
        self.assertFalse(cnr_loess["log2"].isna().any())
        # ...but the actual log2 values differ -- otherwise the smoother
        # choice is being silently ignored somewhere in the call stack.
        self.assertFalse(
            np.allclose(cnr_loess["log2"], cnr_median["log2"]),
            "loess smoother produced identical log2 to median; "
            "the bias_smoother parameter is not reaching center_by_window.",
        )

    def test_fix_bias_smoother_unknown_raises(self):
        """Unknown bias_smoother values raise a clear ValueError, not silently default."""
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        anti_bins = ref[is_bg]
        with self.assertRaises(ValueError) as cm:
            commands.do_fix(tgt_bins, anti_bins, ref, bias_smoother="bogus")
        self.assertIn("bogus", str(cm.exception))

    def test_fix_degenerate_no_nan(self):
        """do_fix never emits NaN log2/weight from degenerate input (#521/#524).

        Zero-coverage sentinel and NaN bins (from malformed inputs or dead
        regions in WGS) must not survive into the .cnr, where they would crash
        smoothing/CBS downstream.
        """
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg].copy()
        anti_bins = ref[is_bg].copy()
        # Inject degenerate values into a few sample bins
        tgt_log2 = tgt_bins["log2"].copy()
        tgt_log2.iloc[0] = np.nan
        tgt_log2.iloc[1] = params.NULL_LOG2_COVERAGE
        tgt_bins["log2"] = tgt_log2
        # Inject a degenerate value into the reference, too
        ref_log2 = ref["log2"].copy()
        ref_log2.iloc[2] = np.nan
        ref["log2"] = ref_log2
        cnr = commands.do_fix(tgt_bins, anti_bins, ref)
        self.assertFalse(np.isnan(cnr["log2"]).any())
        self.assertFalse(np.isnan(cnr["weight"]).any())


class ClusterTests(unittest.TestCase):
    """Tests for clustering used in reference building (cnvlib.cluster)."""

    def test_cluster_hierarchical(self):
        """Test hierarchical clustering on synthetic data with 2 clear groups."""
        from cnvlib.cluster import hierarchical

        rng = np.random.default_rng(42)
        n_bins = 100
        # Group A: 5 samples with similar pattern
        base_a = rng.standard_normal(n_bins)
        group_a = np.array([base_a + rng.normal(0, 0.1, n_bins) for _ in range(5)])
        # Group B: 5 samples with a different pattern
        base_b = rng.standard_normal(n_bins)
        group_b = np.array([base_b + rng.normal(0, 0.1, n_bins) for _ in range(5)])
        samples = np.vstack([group_a, group_b])

        clusters = hierarchical(samples, min_cluster_size=3)
        # Should produce exactly 2 clusters; all samples must be assigned
        self.assertEqual(len(clusters), 2)
        all_indices = sorted(idx for c in clusters for idx in c)
        self.assertEqual(all_indices, list(range(10)))
        # Check they separate the groups correctly
        c0 = set(clusters[0])
        c1 = set(clusters[1])
        group_a_set = set(range(5))
        group_b_set = set(range(5, 10))
        self.assertTrue(
            (c0 == group_a_set and c1 == group_b_set)
            or (c0 == group_b_set and c1 == group_a_set)
        )

    def test_cluster_kmedoids(self):
        """Test k-medoids clustering produces valid clusters."""
        from cnvlib.cluster import kmedoids

        rng = np.random.default_rng(42)
        n_bins = 100
        samples = np.vstack(
            [
                rng.standard_normal((5, n_bins)),
                rng.standard_normal((5, n_bins)) + 3,
            ]
        )
        clusters = kmedoids(samples, k=2)
        self.assertEqual(len(clusters), 2)
        all_indices = sorted(idx for c in clusters for idx in c)
        self.assertEqual(all_indices, list(range(10)))

    def test_create_clusters_column_names(self):
        """Test that create_clusters produces expected column names."""
        rng = np.random.default_rng(42)
        n_bins = 50
        n_samples = 6
        # Create logr_matrix with pseudocount row + sample rows
        logr = np.vstack(
            [
                np.zeros(n_bins),  # pseudocount
                rng.standard_normal((n_samples, n_bins)),
            ]
        )
        depths = np.vstack(
            [
                np.zeros(n_bins),  # pseudocount
                rng.uniform(10, 100, (n_samples, n_bins)),
            ]
        )
        sample_ids = [f"sample_{i}" for i in range(n_samples)]
        cols = reference.create_clusters(
            logr, depths, min_cluster_size=2, sample_ids=sample_ids
        )
        # Should have at least one cluster with log2_, depth_, spread_ columns
        self.assertTrue(len(cols) > 0)
        log2_keys = [k for k in cols if k.startswith("log2_")]
        depth_keys = [k for k in cols if k.startswith("depth_")]
        spread_keys = [k for k in cols if k.startswith("spread_")]
        self.assertTrue(len(log2_keys) > 0)
        self.assertEqual(len(log2_keys), len(depth_keys))
        self.assertEqual(len(log2_keys), len(spread_keys))
        # All column arrays should have n_bins elements
        for val in cols.values():
            self.assertEqual(len(val), n_bins)
