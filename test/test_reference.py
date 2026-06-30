#!/usr/bin/env python
"""Tests for the reference command (and clustering used to build references)."""

import argparse
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


def _write_coverage_cnn(
    path, chrx_log2, chry_log2, chrx_sd, chry_sd, n_auto=200, n_x=40, n_y=40, seed=7
):
    """Write a synthetic coverage .cnn with tunable sex-chromosome signal.

    Used to reproduce target/antitarget sex-inference conflicts: autosomes sit
    at log2 0, and chrX/chrY are drawn around the given means so a caller can
    make, e.g., a haploid chrX with a near-empty (deeply negative) chrY.
    """
    rng = np.random.default_rng(seed)
    rows = [
        (
            chrom,
            i * 1000,
            i * 1000 + 1000,
            "Background",
            100,
            float(rng.normal(0, 0.2)),
        )
        for chrom, n in [("chr1", n_auto), ("chr2", n_auto)]
        for i in range(n)
    ]
    rows.extend(
        (
            "chrX",
            i * 1000,
            i * 1000 + 1000,
            "Background",
            50,
            float(rng.normal(chrx_log2, chrx_sd)),
        )
        for i in range(n_x)
    )
    rows.extend(
        (
            "chrY",
            i * 1000,
            i * 1000 + 1000,
            "Background",
            2,
            float(rng.normal(chry_log2, chry_sd)),
        )
        for i in range(n_y)
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
        """do_fix accepts bias_smoother='loess' end-to-end (#1028).

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

    def test_fix_reference_grammar_argparse(self):
        """`fix` parses the redesigned reference grammar (#894 follow-up).

        The reference moved from a required positional to a canonical
        ``-r/--reference`` flag plus a deprecated optional positional
        (dest ``reference_arg``); the antitarget stays an optional positional.
        This pins how each invocation maps onto the parsed attributes. Crucially,
        the forgotten-reference case (``fix target antitarget``) now parses to an
        *unresolved* reference (both sources None) instead of silently treating
        the antitarget file as the reference -- the silent-misparse bug this
        redesign closes. argparse no longer enforces the reference (that moved
        into ``_cmd_fix``), so none of these forms raise SystemExit.
        """
        # WGS form: target + -r reference, antitarget omitted.
        ns = commands.P_fix.parse_args(["tgt.cnn", "-r", "ref.cnn"])
        self.assertEqual(ns.target, "tgt.cnn")
        self.assertIsNone(ns.antitarget)
        self.assertIsNone(ns.reference_arg)
        self.assertEqual(ns.reference, "ref.cnn")
        # Hybrid-capture form: target + antitarget positional + -r reference.
        ns = commands.P_fix.parse_args(["tgt.cnn", "anti.cnn", "-r", "ref.cnn"])
        self.assertEqual(ns.target, "tgt.cnn")
        self.assertEqual(ns.antitarget, "anti.cnn")
        self.assertIsNone(ns.reference_arg)
        self.assertEqual(ns.reference, "ref.cnn")
        # Deprecated 3-positional form: reference lands in reference_arg, not -r.
        ns = commands.P_fix.parse_args(["tgt.cnn", "anti.cnn", "ref.cnn"])
        self.assertEqual(ns.target, "tgt.cnn")
        self.assertEqual(ns.antitarget, "anti.cnn")
        self.assertEqual(ns.reference_arg, "ref.cnn")
        self.assertIsNone(ns.reference)
        # Forgotten-reference trap: the old grammar took this as target+reference;
        # the new one parses it as target+antitarget with the reference left
        # unresolved (both reference sources None), to be rejected by _cmd_fix.
        ns = commands.P_fix.parse_args(["tgt.cnn", "anti.cnn"])
        self.assertEqual(ns.target, "tgt.cnn")
        self.assertEqual(ns.antitarget, "anti.cnn")
        self.assertIsNone(ns.reference_arg)
        self.assertIsNone(ns.reference)
        # Lone target: antitarget and both reference sources unset.
        ns = commands.P_fix.parse_args(["tgt.cnn"])
        self.assertEqual(ns.target, "tgt.cnn")
        self.assertIsNone(ns.antitarget)
        self.assertIsNone(ns.reference_arg)
        self.assertIsNone(ns.reference)

    @staticmethod
    def _fix_namespace(
        target,
        *,
        antitarget_fname=None,
        reference=None,
        reference_arg=None,
        output=None,
    ):
        """Build a `fix` argparse.Namespace, varying only the reference inputs.

        do_edge is held False so callers need not assert value equality against
        do_fix's do_edge=True default; only bin counts and validity are compared.
        """
        return argparse.Namespace(
            target=target,
            antitarget=antitarget_fname,
            reference_arg=reference_arg,
            reference=reference,
            sample_id=None,
            diploid_parx_genome=None,
            do_gc=True,
            do_edge=False,
            do_rmask=True,
            cluster=False,
            smoothing_window_fraction=None,
            bias_smoother="median",
            output=output,
        )

    def test_cmd_fix_without_antitarget(self):
        """`_cmd_fix` produces a valid .cnr from the WGS `-r` form (#894).

        Drives the full command path a WGS/amplicon user takes (target plus
        ``-r reference``, no antitarget file): reference via the canonical
        ``reference`` attribute, ``reference_arg`` None, ``antitarget`` None.
        The None-antitarget branch must synthesize an empty antitarget and still
        write a usable .cnr -- not crash and not emit NaN log2 (which would break
        downstream segmentation). Cross-checks that this yields the same number
        of bins as the explicit blank-antitarget API call, confirming the
        synthesized empty antitarget is equivalent to passing one.
        """
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        tmpdir = tempfile.mkdtemp()
        try:
            tgt_fname = os.path.join(tmpdir, "sample.targetcoverage.cnn")
            ref_fname = os.path.join(tmpdir, "reference.cnn")
            out_fname = os.path.join(tmpdir, "sample.cnr")
            tabio.write(tgt_bins, tgt_fname)
            tabio.write(ref[~is_bg], ref_fname)
            namespace = self._fix_namespace(
                tgt_fname, reference=ref_fname, output=out_fname
            )
            commands._cmd_fix(namespace)
            cnr = cnvlib.read(out_fname)
            # (a) Produced real bins, bounded by the target input.
            self.assertTrue(0 < len(cnr) <= len(tgt_bins))
            # (b) No NaN log2 leaked through (would crash downstream segmentation).
            self.assertFalse(np.isnan(cnr["log2"]).any())
            # (c) Same bin count as the explicit blank-antitarget API path.
            #     (Value equality is intentionally not asserted: this namespace
            #     sets do_edge=False while do_fix defaults to do_edge=True, so
            #     the corrected log2 values legitimately differ.)
            api_cnr = commands.do_fix(
                tgt_bins, cnvlib.cnary.CopyNumArray([]), ref[~is_bg]
            )
            self.assertEqual(len(cnr), len(api_cnr))
        finally:
            shutil.rmtree(tmpdir)

    def test_cmd_fix_requires_reference(self):
        """`_cmd_fix` exits when no reference is given by flag or positional (#894).

        The forgotten-reference invocation (`fix target antitarget`) parses with
        both reference sources None; _cmd_fix must refuse rather than proceed
        referenceless. A real target file is provided so the exit is the
        reference check, not file I/O.
        """
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        tmpdir = tempfile.mkdtemp()
        try:
            tgt_fname = os.path.join(tmpdir, "sample.targetcoverage.cnn")
            tabio.write(tgt_bins, tgt_fname)
            namespace = self._fix_namespace(
                tgt_fname, reference=None, reference_arg=None
            )
            with self.assertRaisesRegex(SystemExit, "No reference"):
                commands._cmd_fix(namespace)
        finally:
            shutil.rmtree(tmpdir)

    def test_cmd_fix_reference_specified_twice(self):
        """`_cmd_fix` exits when the reference is given both ways (#894).

        Supplying -r/--reference *and* the deprecated positional is ambiguous;
        the command must refuse rather than silently pick one.
        """
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        tmpdir = tempfile.mkdtemp()
        try:
            tgt_fname = os.path.join(tmpdir, "sample.targetcoverage.cnn")
            ref_fname = os.path.join(tmpdir, "reference.cnn")
            tabio.write(tgt_bins, tgt_fname)
            tabio.write(ref[~is_bg], ref_fname)
            namespace = self._fix_namespace(
                tgt_fname, reference=ref_fname, reference_arg=ref_fname
            )
            with self.assertRaisesRegex(SystemExit, "once"):
                commands._cmd_fix(namespace)
        finally:
            shutil.rmtree(tmpdir)

    def test_cmd_fix_positional_reference_deprecated(self):
        """The legacy positional reference still works but warns (#894).

        Passing the reference as the third positional (reference_arg) must emit a
        DeprecationWarning yet remain fully functional -- the deprecation is
        non-fatal, so a valid .cnr is still written. Uses the full hybrid fixture
        (target + antitarget + full reference) so do_fix succeeds.
        """
        ref = cnvlib.read("formats/reference-tr.cnn")
        is_bg = ref["gene"] == "Background"
        tgt_bins = ref[~is_bg]
        anti_bins = ref[is_bg]
        tmpdir = tempfile.mkdtemp()
        try:
            tgt_fname = os.path.join(tmpdir, "sample.targetcoverage.cnn")
            anti_fname = os.path.join(tmpdir, "sample.antitargetcoverage.cnn")
            ref_fname = os.path.join(tmpdir, "reference.cnn")
            out_fname = os.path.join(tmpdir, "sample.cnr")
            tabio.write(tgt_bins, tgt_fname)
            tabio.write(anti_bins, anti_fname)
            tabio.write(ref, ref_fname)
            namespace = self._fix_namespace(
                tgt_fname,
                antitarget_fname=anti_fname,
                reference_arg=ref_fname,
                output=out_fname,
            )
            with self.assertWarns(DeprecationWarning):
                commands._cmd_fix(namespace)
            cnr = cnvlib.read(out_fname)
            self.assertTrue(0 < len(cnr) <= len(ref))
            self.assertFalse(np.isnan(cnr["log2"]).any())
        finally:
            shutil.rmtree(tmpdir)

    def test_match_ref_to_sample_zero_overlap(self):
        """Zero shared bins reads as a wrong-file error, not a panel mismatch.

        When the reference and sample share *no* coordinates (e.g. a raw coverage
        .cnn handed in where a built reference was expected, or the reference
        omitted so the antitarget file gets misused), match_ref_to_sample must
        say so explicitly rather than emit a misleading "missing N bins".
        """
        samp = cnary.CopyNumArray.from_columns(
            {
                "chromosome": ["chr1"] * 4,
                "start": [0, 100, 200, 300],
                "end": [100, 200, 300, 400],
                "gene": list("ABCD"),
                "log2": [0.0, 0.1, -0.1, 0.2],
            },
            {"sample_id": "Sample"},
        )
        # Disjoint coordinates: a different chromosome -> no bins in common.
        ref = cnary.CopyNumArray.from_columns(
            {
                "chromosome": ["chr2"] * 4,
                "start": [0, 100, 200, 300],
                "end": [100, 200, 300, 400],
                "gene": list("ABCD"),
                "log2": [0.0, 0.0, 0.0, 0.0],
            }
        )
        with self.assertRaisesRegex(ValueError, "shares no bins"):
            fix.match_ref_to_sample(ref, samp)

    def test_match_ref_to_sample_partial_overlap(self):
        """A reference missing *some* of the sample's bins flags a panel mismatch.

        Partial overlap (the reference covers part but not all of the sample's
        bins) is the "wrong target panel" case and must be distinguished from the
        zero-overlap wrong-file case: the message reports how many of how many
        bins are missing.
        """
        samp = cnary.CopyNumArray.from_columns(
            {
                "chromosome": ["chr1"] * 4,
                "start": [0, 100, 200, 300],
                "end": [100, 200, 300, 400],
                "gene": list("ABCD"),
                "log2": [0.0, 0.1, -0.1, 0.2],
            },
            {"sample_id": "Sample"},
        )
        # Reference covers only the first two of the sample's four bins.
        ref = cnary.CopyNumArray.from_columns(
            {
                "chromosome": ["chr1", "chr1"],
                "start": [0, 100],
                "end": [100, 200],
                "gene": ["A", "B"],
                "log2": [0.0, 0.0],
            }
        )
        with self.assertRaisesRegex(ValueError, r"missing \d+ of \d+"):
            fix.match_ref_to_sample(ref, samp)


class ClusterTests(unittest.TestCase):
    """Tests for clustering used in reference building (cnvlib.cluster)."""

    def test_cluster_hierarchical(self):
        """Test hierarchical clustering on synthetic data with 2 clear groups."""
        rng = np.random.default_rng(42)
        n_bins = 100
        # Group A: 5 samples with similar pattern
        base_a = rng.standard_normal(n_bins)
        group_a = np.array([base_a + rng.normal(0, 0.1, n_bins) for _ in range(5)])
        # Group B: 5 samples with a different pattern
        base_b = rng.standard_normal(n_bins)
        group_b = np.array([base_b + rng.normal(0, 0.1, n_bins) for _ in range(5)])
        samples = np.vstack([group_a, group_b])

        clusters = cluster.hierarchical(samples, min_cluster_size=3)
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
        rng = np.random.default_rng(42)
        n_bins = 100
        samples = np.vstack(
            [
                rng.standard_normal((5, n_bins)),
                rng.standard_normal((5, n_bins)) + 3,
            ]
        )
        clusters = cluster.kmedoids(samples, k=2)
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
