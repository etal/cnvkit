#!/usr/bin/env python
"""Tests for the batch command."""

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


class BatchTests(unittest.TestCase):
    """Tests for the batch command and end-to-end pipeline."""

    def test_batch(self):
        """The 'batch' command."""
        target_bed = "formats/my-targets.bed"
        fasta = "formats/chrM-Y-trunc.hg19.fa"
        bam = "formats/na12878-chrM-Y-trunc.bam"
        annot = "formats/my-refflat.bed"
        # Build a single-sample WGS reference
        ref_fname, tgt_bed_fname, _ = batch.batch_make_reference(
            [bam],
            None,
            None,
            True,
            None,
            fasta,
            annot,
            True,
            500,
            None,
            None,
            None,
            None,
            "build",
            1,
            False,
            0,
            "wgs",
            False,
        )
        self.assertEqual(ref_fname, "build/reference.cnn")
        refarr = cnvlib.read(ref_fname, "bed")
        tgt_regions = tabio.read(tgt_bed_fname, "bed")
        self.assertEqual(len(refarr), len(tgt_regions))
        # Build a single-sample hybrid-capture reference
        ref_fname, tgt_bed_fname, anti_bed_fname = batch.batch_make_reference(
            [bam],
            target_bed,
            None,
            True,
            None,
            fasta,
            None,
            True,
            10,
            None,
            1000,
            100,
            None,
            "build",
            1,
            False,
            0,
            "hybrid",
            False,
        )
        self.assertEqual(ref_fname, "build/reference.cnn")
        refarr = cnvlib.read(ref_fname, "bed")
        tgt_regions = tabio.read(tgt_bed_fname, "bed")
        anti_regions = tabio.read(anti_bed_fname, "bed")
        self.assertEqual(len(refarr), len(tgt_regions) + len(anti_regions))
        # Run the same sample
        batch.batch_run_sample(
            bam,
            tgt_bed_fname,
            anti_bed_fname,
            ref_fname,
            "build",
            True,
            None,
            True,
            True,
            "Rscript",
            False,
            0,
            False,
            "hybrid",
            "hmm",
            1,
            False,
        )
        cns = cnvlib.read("build/na12878-chrM-Y-trunc.cns")
        self.assertGreater(len(cns), 0)

    def test_batch_hybrid_autobin(self):
        """Hybrid batch with autobin-determined bin sizes (issue #302)."""
        target_bed = "formats/my-targets.bed"
        fasta = "formats/chrM-Y-trunc.hg19.fa"
        bam = "formats/na12878-chrM-Y-trunc.bam"
        # target_avg_size=0 triggers autobin for hybrid mode
        ref_fname, tgt_bed_fname, anti_bed_fname = batch.batch_make_reference(
            [bam],
            target_bed,
            None,
            True,
            None,
            fasta,
            None,
            True,
            0,
            None,
            None,
            None,
            None,
            "build",
            1,
            False,
            0,
            "hybrid",
            False,
        )
        self.assertEqual(ref_fname, "build/reference.cnn")
        refarr = cnvlib.read(ref_fname, "bed")
        tgt_regions = tabio.read(tgt_bed_fname, "bed")
        anti_regions = tabio.read(anti_bed_fname, "bed")
        self.assertGreater(len(tgt_regions), 0)
        self.assertEqual(len(refarr), len(tgt_regions) + len(anti_regions))

    def test_batch_run_sample_with_duplicate_coordinates(self):
        """Integration test: batch_run_sample with duplicate coordinates in coverage.

        This tests that when coverage data contains duplicate coordinates,
        the error is properly detected and reported with sample context.
        """
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create mock coverage data with duplicate coordinates
            # This simulates the issue from GitHub issue #971
            target_data = pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 4,
                    "start": [100, 200, 200, 300],  # Duplicate start at position 200
                    "end": [150, 250, 250, 350],  # Duplicate end at position 250
                    "gene": ["GeneA", "GeneB", "GeneB", "GeneC"],
                    "log2": [0.1, 0.2, 0.2, 0.3],
                    "depth": [100.0, 200.0, 200.0, 150.0],
                }
            )
            target_cnarr = cnary.CopyNumArray(target_data, {"sample_id": "test_sample"})

            antitarget_data = pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 2,
                    "start": [50, 400],
                    "end": [90, 450],
                    "gene": ["Antitarget", "Antitarget"],
                    "log2": [0.0, 0.0],
                    "depth": [50.0, 50.0],
                }
            )
            antitarget_cnarr = cnary.CopyNumArray(
                antitarget_data, {"sample_id": "test_sample"}
            )

            # Create a simple reference without duplicates
            ref_data = pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 5,
                    "start": [50, 100, 200, 300, 400],
                    "end": [90, 150, 250, 350, 450],
                    "gene": ["Antitarget", "GeneA", "GeneB", "GeneC", "Antitarget"],
                    "log2": [0.0, 0.0, 0.0, 0.0, 0.0],
                }
            )
            ref_cnarr = cnary.CopyNumArray(ref_data, {"sample_id": "reference"})

            # Test that do_fix detects duplicate coordinates and raises ValueError
            with self.assertRaises(ValueError) as ctx:
                fix.do_fix(target_cnarr, antitarget_cnarr, ref_cnarr)

            # The error message should mention duplicates
            self.assertIn("Duplicated genomic coordinates", str(ctx.exception))

    def test_batch_accepts_sample_sex_arg(self):
        """``cnvkit.py batch`` accepts ``-x/--sample-sex`` (#500, #635).

        Users running the full pipeline via ``batch`` need a way to override
        the auto-inferred sex when the inference goes wrong on their data,
        without having to hand-build the equivalent multi-step pipeline.
        The ``call``/``diagram``/``export`` commands have always accepted
        ``--sample-sex``; ``batch`` was the documented gap.
        """
        args = commands.AP.parse_args(
            [
                "batch",
                "dummy.bam",
                "-r",
                "dummy.cnn",
                "-d",
                "outdir",
                "-x",
                "female",
            ]
        )
        self.assertEqual(args.sample_sex, "female")
        # Long form also resolves to sample_sex. (The legacy ``-g/--gender``
        # alias supplied by ``add_sample_sex`` elsewhere is deliberately
        # omitted on batch because ``-g`` is already taken by ``--access``.)
        args2 = commands.AP.parse_args(
            [
                "batch",
                "dummy.bam",
                "-r",
                "dummy.cnn",
                "-d",
                "outdir",
                "--sample-sex",
                "male",
            ]
        )
        self.assertEqual(args2.sample_sex, "male")

    def test_batch_run_sample_propagates_sex_to_do_call(self):
        """``batch_run_sample`` must pass ``is_haploid_x_reference`` AND
        ``is_sample_female`` to every ``call.do_call`` invocation, not silently
        fall through to ``do_call``'s ``False/False`` defaults.

        Previously the ``do_call`` calls inside ``batch_run_sample`` used the
        function defaults regardless of ``batch -y``, so ``.call.cns`` chrX
        cn values were computed against a diploid-X reference even when the
        user built a male reference (normal female chrX at log2 ~ +1 was
        called as GAIN). The fix consolidates the per-sample sex decision
        via ``verify_sample_sex`` and passes the result explicitly to both
        ``do_call`` calls and the diagram path.

        AST-level guard so source rearrangements still catch the regression
        cleanly without needing a full BAM fixture per case.
        """
        import ast
        import inspect

        src = inspect.getsource(batch.batch_run_sample)
        tree = ast.parse(src)
        do_call_invocations = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "do_call"
                and isinstance(func.value, ast.Name)
                and func.value.id == "call"
            ):
                do_call_invocations.append(node)
        self.assertGreater(
            len(do_call_invocations),
            0,
            "Expected at least one call.do_call() invocation in batch_run_sample",
        )
        for inv in do_call_invocations:
            kw_names = {kw.arg for kw in inv.keywords}
            self.assertIn(
                "is_haploid_x_reference",
                kw_names,
                "call.do_call inside batch_run_sample must pass "
                "is_haploid_x_reference explicitly (#500/#635-adjacent latent bug).",
            )
            self.assertIn(
                "is_sample_female",
                kw_names,
                "call.do_call inside batch_run_sample must pass "
                "is_sample_female explicitly so batch -x is honored.",
            )

    def test_batch_run_sample_passes_bias_smoother_to_do_fix(self):
        """``batch_run_sample`` must forward its ``bias_smoother`` argument
        to every ``fix.do_fix`` invocation it makes (gh#1028).

        Without this, ``cnvkit batch --bias-smoother loess`` would silently
        fall back to ``rolling_median`` and users would have no way to
        evaluate LOESS through the high-level ``batch`` entry point even
        though the CLI flag was accepted.

        AST-level guard so source rearrangements still catch the regression
        cleanly without needing a full BAM fixture.
        """
        import ast
        import inspect

        src = inspect.getsource(batch.batch_run_sample)
        tree = ast.parse(src)
        do_fix_invocations = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "do_fix"
                and isinstance(func.value, ast.Name)
                and func.value.id == "fix"
            ):
                do_fix_invocations.append(node)
        self.assertGreater(
            len(do_fix_invocations),
            0,
            "Expected at least one fix.do_fix() invocation in batch_run_sample",
        )
        for inv in do_fix_invocations:
            kw_names = {kw.arg for kw in inv.keywords}
            self.assertIn(
                "bias_smoother",
                kw_names,
                "fix.do_fix inside batch_run_sample must pass bias_smoother "
                "explicitly so `cnvkit batch --bias-smoother loess` reaches "
                "center_by_window (gh#1028).",
            )

    def test_batch_make_reference_passes_bias_smoother_to_do_reference(self):
        """``batch_make_reference`` must forward ``bias_smoother`` to
        ``reference.do_reference`` for the pooled-reference path (gh#1028).

        ``do_reference_flat`` is exempt because a flat reference performs no
        bias-vs-trait smoothing.

        AST-level guard, paired with
        ``test_batch_run_sample_passes_bias_smoother_to_do_fix``.
        """
        import ast
        import inspect

        src = inspect.getsource(batch.batch_make_reference)
        tree = ast.parse(src)
        do_reference_invocations = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "do_reference"
                and isinstance(func.value, ast.Name)
                and func.value.id == "reference"
            ):
                do_reference_invocations.append(node)
        self.assertGreater(
            len(do_reference_invocations),
            0,
            "Expected at least one reference.do_reference() invocation "
            "in batch_make_reference",
        )
        for inv in do_reference_invocations:
            kw_names = {kw.arg for kw in inv.keywords}
            self.assertIn(
                "bias_smoother",
                kw_names,
                "reference.do_reference inside batch_make_reference must "
                "pass bias_smoother explicitly so `cnvkit batch "
                "--bias-smoother loess` reaches center_by_window during "
                "reference construction (gh#1028).",
            )

    @pytest.mark.slow
    def test_diploid_parx_genome(self):
        genome_build = "grch37"
        male_target_fnames = [
            "formats/male01.targetcoverage.cnn",
            "formats/male02.targetcoverage.cnn",
        ]
        male_antitarget_fnames = [
            "formats/male01.antitargetcoverage.cnn",
            "formats/male02.antitargetcoverage.cnn",
        ]
        female_target_fnames = [
            "formats/female01.targetcoverage.cnn",
            "formats/female02.targetcoverage.cnn",
        ]
        female_antitarget_fnames = [
            "formats/female01.antitargetcoverage.cnn",
            "formats/female02.antitargetcoverage.cnn",
        ]

        # Baseline assumption: PAR1/2 is highly coveraged in the male sample.
        target_cnn = cnvlib.read(male_target_fnames[0])
        antitarget_cnn = cnvlib.read(male_antitarget_fnames[0])
        target_cnn_non_parx = target_cnn[target_cnn.chr_x_filter(genome_build)]
        antitarget_cnn_non_parx = antitarget_cnn[
            antitarget_cnn.chr_x_filter(genome_build)
        ]
        target_cnn_parx = target_cnn[target_cnn.parx_filter(genome_build)]
        antitarget_cnn_parx = antitarget_cnn[antitarget_cnn.parx_filter(genome_build)]
        target_log2_mean_parx = target_cnn_parx["log2"].mean()
        target_log2_mean_non_parx = target_cnn_non_parx["log2"].mean()
        antitarget_log2_mean_parx = antitarget_cnn_parx["log2"].mean()
        antitarget_log2_mean_non_parx = antitarget_cnn_non_parx["log2"].mean()
        self.assertTrue(
            target_log2_mean_parx - 0.9 > target_log2_mean_non_parx,
            "PAR1/2 have nearly doubled coverage.",
        )
        self.assertTrue(
            antitarget_log2_mean_parx - 0.9 > antitarget_log2_mean_non_parx,
            "PAR1/2 have nearly doubled coverage.",
        )

        #### MIXED POOL ####
        target_fnames = male_target_fnames + female_target_fnames
        antitarget_fnames = male_antitarget_fnames + female_antitarget_fnames
        # Only process male02 (index 1) and female01 (index 2) — the samples
        # actually checked below — to avoid redundant CBS segmentation.
        _ref_probes1, _cnrs1, _cnss1, clls1, _sex_df1 = run_samples(
            target_fnames, antitarget_fnames, None, sample_indices=[1, 2]
        )
        _ref_probes2, _cnrs2, _cnss2, clls2, _sex_df2 = run_samples(
            target_fnames, antitarget_fnames, genome_build, sample_indices=[1, 2]
        )

        # "dpxg" = DiploidParXGenome
        male_call, male_call_dpxg = clls1[0], clls2[0]
        female_call, female_call_dpxg = clls1[1], clls2[1]
        male_call__x, male_call_dpxg__x = (
            male_call[male_call.chr_x_filter()],
            male_call_dpxg[male_call_dpxg.chromosome == "X"],
        )
        female_call__x, female_call_dpxg__x = (
            female_call[female_call.chromosome == "X"],
            female_call_dpxg[female_call_dpxg.chromosome == "X"],
        )

        # The cn=0 segment for the male sample is derived from a true loss on XAGE1B.
        self.assertEqual(
            male_call__x["cn"].to_list(),
            [1, 1, 0, 1, 1],
            "Non-Dpxg male has cn=1 including PAR1/2.",
        )
        self.assertEqual(
            male_call_dpxg__x["cn"].to_list(),
            [2, 1, 0, 1, 1, 2],
            "Dpxg male has cn=2 for PAR1/2 and cn=1 otherwise.",
        )
        self.assertEqual(
            female_call__x["cn"].to_list(),
            [1, 2, 2, 1],
            "Non-Dpxg female is biased towards loss in PAR1/2.",
        )
        self.assertEqual(
            female_call_dpxg__x["cn"].to_list(),
            [2, 2],
            "Dpxg female has cn=2 everywhere including PAR1/2.",
        )

        #### MALE-ONLY POOL ####
        target_fnames = male_target_fnames
        antitarget_fnames = male_antitarget_fnames
        # Only process male02 (index 1) — the sample actually checked below.
        _ref_probes1, _cnrs1, _cnss1, clls1, _sex_df1 = run_samples(
            target_fnames, antitarget_fnames, None, sample_indices=[1]
        )
        _ref_probes2, _cnrs2, _cnss2, clls2, _sex_df2 = run_samples(
            target_fnames, antitarget_fnames, genome_build, sample_indices=[1]
        )

        male_call, male_call_dpxg = clls1[0], clls2[0]
        male_call__x, male_call_dpxg__x = (
            male_call[male_call.chr_x_filter()],
            male_call_dpxg[male_call_dpxg.chromosome == "X"],
        )
        self.assertEqual(
            male_call__x["cn"].to_list(),
            [1, 1],
            "Non-Dpxg male has cn=1 including PAR1/2.",
        )
        self.assertEqual(
            male_call_dpxg__x["cn"].to_list(),
            [2, 1, 1, 2],
            "Dpxg male has cn=2 for PAR1/2 and cn=1 otherwise.",
        )


# == helpers ==


def run_samples(
    target_fnames, antitarget_fnames, diploid_parx_genome, sample_indices=None
):
    ref_probes = commands.do_reference(
        target_fnames, antitarget_fnames, diploid_parx_genome=diploid_parx_genome
    )
    if sample_indices is None:
        sample_indices = range(len(target_fnames))
    cnrs, cnss, clls = [], [], []
    for i in sample_indices:
        cnr, cns, cll = run_sample(
            target_fnames[i], antitarget_fnames[i], ref_probes, diploid_parx_genome
        )
        cnrs.append(cnr)
        cnss.append(cns)
        clls.append(cll)

    sex_df = commands.do_sex(
        cnrs, is_haploid_x_reference=False, diploid_parx_genome=diploid_parx_genome
    )
    return ref_probes, cnrs, cnss, clls, sex_df


def run_sample(target_fname, antitarget_fname, ref_probes, diploid_parx_genome):
    tgt_raw = cnvlib.read(target_fname)
    anti_raw = cnvlib.read(antitarget_fname)
    cnr = commands.do_fix(
        tgt_raw, anti_raw, ref_probes, diploid_parx_genome=diploid_parx_genome
    )
    cns = commands.do_segmentation(
        cnr, method="cbs", diploid_parx_genome=diploid_parx_genome, threshold=0.001
    )
    cll = commands.do_call(cns, diploid_parx_genome=diploid_parx_genome)
    return (cnr, cns, cll)


if __name__ == "__main__":
    unittest.main(verbosity=2)
