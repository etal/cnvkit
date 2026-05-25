#!/usr/bin/env python
"""Tests for the call command and heterozygous-SNP loading."""

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


class CallTests(unittest.TestCase):
    """Tests for the call command and related functions."""

    def test_call(self):
        """The 'call' command."""
        # Methods: clonal, threshold, none
        tr_cns = cnvlib.read("formats/tr95t.cns")
        tr_thresh = commands.do_call(
            tr_cns,
            None,
            "threshold",
            is_haploid_x_reference=True,
            is_sample_female=True,
        )
        self.assertEqual(len(tr_cns), len(tr_thresh))
        tr_clonal = commands.do_call(
            tr_cns,
            None,
            "clonal",
            purity=0.65,
            is_haploid_x_reference=True,
            is_sample_female=True,
        )
        self.assertEqual(len(tr_cns), len(tr_clonal))
        cl_cns = cnvlib.read("formats/cl_seq.cns")
        cl_thresh = commands.do_call(
            cl_cns,
            None,
            "threshold",
            thresholds=np.log2((np.arange(12) + 0.5) / 6.0),
            is_haploid_x_reference=True,
            is_sample_female=True,
        )
        self.assertEqual(len(cl_cns), len(cl_thresh))
        cl_clonal = commands.do_call(
            cl_cns,
            None,
            "clonal",
            ploidy=6,
            purity=0.99,
            is_haploid_x_reference=True,
            is_sample_female=True,
        )
        self.assertEqual(len(cl_cns), len(cl_clonal))
        cl_none = commands.do_call(
            cl_cns,
            None,
            "none",
            ploidy=6,
            purity=0.99,
            is_haploid_x_reference=True,
            is_sample_female=True,
        )
        self.assertEqual(len(cl_cns), len(cl_none))

    def test_call_filter(self):
        segments = cnvlib.read("formats/tr95t.segmetrics.cns")
        variants = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        # Each filter individually, then filter combos (with merges where needed)
        for filters, merges in (
            (["ampdel"], []),
            (["cn"], []),
            (["ci"], []),
            (["sem"], []),
            ([], ["bic"]),
            (["sem", "cn", "ampdel"], []),
            (["ci", "cn"], []),
            (["ci", "sem"], []),
            (["cn", "ci", "sem"], []),
            ([], ["bic", "cn"]),
            (["ci"], ["bic"]),
            (["ci", "cn"], ["bic"]),
        ):
            with self.subTest(filters=filters, merges=merges):
                result = commands.do_call(
                    segments,
                    variants,
                    method="threshold",
                    purity=0.9,
                    is_haploid_x_reference=True,
                    is_sample_female=True,
                    filters=filters or None,
                    merges=merges or None,
                )
                self.assertLessEqual(len(result), len(segments))
                if "ampdel" not in filters:
                    # At least 1 segment per chromosome remains
                    self.assertLessEqual(len(segments.chromosome.unique()), len(result))
                for colname in "baf", "cn", "cn1", "cn2":
                    self.assertIn(colname, result)
                # Segmetrics columns must survive merging by prior filters
                for colname in "sem", "ci_lo", "ci_hi", "stdev":
                    self.assertIn(
                        colname,
                        result,
                        f"{colname!r} dropped after filters={filters} merges={merges}",
                    )

    def test_call_filter_ci_preserves_different_magnitudes(self):
        """CI filter should not merge adjacent segments with different magnitudes.

        Two adjacent segments both with CIs entirely below zero (both are
        confident losses) should remain separate if they have different log2
        values.  Only segments whose CI overlaps zero (neutral) should be
        merged with adjacent neutral segments.
        """
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 5,
                    "start": [0, 1000, 2000, 3000, 4000],
                    "end": [1000, 2000, 3000, 4000, 5000],
                    "gene": ["A", "B", "C", "D", "E"],
                    "log2": [-1.0, -0.03, 0.01, 0.02, 0.5],
                    "weight": [1.0, 1.0, 1.0, 1.0, 1.0],
                    "probes": [10, 10, 10, 10, 10],
                    # Segment A: deep loss, CI entirely below 0
                    # Segment B: slight loss, CI entirely below 0
                    # Segment C: neutral, CI overlaps 0
                    # Segment D: neutral, CI overlaps 0
                    # Segment E: gain, CI entirely above 0
                    "ci_lo": [-1.2, -0.06, -0.05, -0.04, 0.3],
                    "ci_hi": [-0.8, -0.01, 0.07, 0.08, 0.7],
                }
            )
        )

        result = segfilters.ci(segarr)

        # A and B must stay separate -- both are losses but different magnitude
        # C and D should merge -- both are neutral (CI overlaps zero)
        # E must stay separate -- it's a gain
        self.assertEqual(len(result), 4)
        # Check that the deep loss is preserved with its own log2
        self.assertAlmostEqual(result["log2"].iat[0], -1.0)
        # Check that the slight loss is preserved with its own log2
        self.assertAlmostEqual(result["log2"].iat[1], -0.03)
        # Check that the two neutrals were merged
        self.assertEqual(result["probes"].iat[2], 20)
        # Check that the gain is preserved
        self.assertAlmostEqual(result["log2"].iat[3], 0.5)

    def test_call_filter_bic_different_means_stay_separate(self):
        """BIC filter should not merge segments with clearly different means."""
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 3,
                    "start": [0, 1000, 2000],
                    "end": [1000, 2000, 3000],
                    "gene": ["A", "B", "C"],
                    "log2": [-1.0, 0.0, 0.5],
                    "weight": [10.0, 10.0, 10.0],
                    "probes": [50, 50, 50],
                    "stdev": [0.1, 0.1, 0.1],
                }
            )
        )
        result = segfilters.bic(segarr)
        # All three segments have very different means with low variance
        self.assertEqual(len(result), 3)

    def test_call_filter_bic_similar_means_merge(self):
        """BIC filter should merge adjacent segments with similar means."""
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 3,
                    "start": [0, 1000, 2000],
                    "end": [1000, 2000, 3000],
                    "gene": ["A", "B", "C"],
                    "log2": [0.31, 0.30, 0.29],
                    "weight": [10.0, 10.0, 10.0],
                    "probes": [20, 20, 20],
                    "stdev": [0.5, 0.5, 0.5],
                }
            )
        )
        result = segfilters.bic(segarr)
        # All three are essentially the same; BIC should merge them
        self.assertEqual(len(result), 1)
        self.assertEqual(result["probes"].iat[0], 60)

    def test_call_filter_bic_chromosome_boundary(self):
        """BIC filter should never merge segments on different chromosomes."""
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1", "chr2"],
                    "start": [0, 0],
                    "end": [1000, 1000],
                    "gene": ["A", "B"],
                    "log2": [0.3, 0.3],
                    "weight": [10.0, 10.0],
                    "probes": [20, 20],
                    "stdev": [0.5, 0.5],
                }
            )
        )
        result = segfilters.bic(segarr)
        # Same means but different chromosomes: must not merge
        self.assertEqual(len(result), 2)

    def test_call_filter_bic_single_bin_fallback(self):
        """BIC filter handles segments with stdev=0 via fallback variance."""
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 3,
                    "start": [0, 1000, 2000],
                    "end": [1000, 2000, 3000],
                    "gene": ["A", "B", "C"],
                    "log2": [0.3, 0.3, 0.3],
                    "weight": [1.0, 10.0, 10.0],
                    "probes": [1, 20, 20],
                    "stdev": [0.0, 0.5, 0.5],
                }
            )
        )
        result = segfilters.bic(segarr)
        # Segment A has stdev=0 (single bin), fallback variance used;
        # all have similar means so they should merge
        self.assertEqual(len(result), 1)
        self.assertEqual(result["probes"].iat[0], 41)

    def test_call_filter_bic_cascading_merges(self):
        """BIC filter iteratively merges: A~B and B~C yields A+B+C."""
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 4,
                    "start": [0, 1000, 2000, 3000],
                    "end": [1000, 2000, 3000, 4000],
                    "gene": ["A", "B", "C", "D"],
                    "log2": [0.10, 0.12, 0.14, -0.5],
                    "weight": [10.0, 10.0, 10.0, 10.0],
                    "probes": [20, 20, 20, 20],
                    "stdev": [0.4, 0.4, 0.4, 0.1],
                }
            )
        )
        result = segfilters.bic(segarr)
        # A, B, C have similar means and high variance -> merge into one
        # D has a clearly different mean -> stays separate
        self.assertEqual(len(result), 2)
        self.assertEqual(result["probes"].iat[0], 60)
        self.assertEqual(result["probes"].iat[1], 20)

    def test_call_merge_bic_via_do_call(self):
        """BIC merge works through the do_call pipeline."""
        segments = cnvlib.read("formats/tr95t.segmetrics.cns")
        result = commands.do_call(
            segments,
            variants=None,
            method="threshold",
            purity=0.9,
            is_haploid_x_reference=True,
            is_sample_female=True,
            merges=["bic"],
        )
        self.assertGreater(len(result), 0)
        self.assertLessEqual(len(result), len(segments))

    def test_call_filter_cn_vs_merge_cn(self):
        """--filter cn merges only neutral segments; --merge cn merges all same-CN."""
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 4,
                    "start": [0, 1000, 2000, 3000],
                    "end": [1000, 2000, 3000, 4000],
                    "gene": ["A", "B", "C", "D"],
                    "log2": [-0.5, -0.6, 0.0, 0.01],
                    "weight": [1.0, 1.0, 1.0, 1.0],
                    "probes": [10, 10, 10, 10],
                }
            )
        )
        # --filter cn: only merges adjacent neutral (cn==2) segments
        result_filter = commands.do_call(
            segarr,
            method="threshold",
            filters=["cn"],
        )
        # Segments A,B are cn=1 (losses) — NOT merged by --filter cn
        # Segments C,D are cn=2 (neutral) — merged by --filter cn
        self.assertEqual(len(result_filter), 3)

        # --merge cn: merges any adjacent same-CN segments
        result_merge = commands.do_call(
            segarr,
            method="threshold",
            merges=["cn"],
        )
        # Segments A,B are cn=1 — merged by --merge cn
        # Segments C,D are cn=2 — merged by --merge cn
        self.assertEqual(len(result_merge), 2)

    def test_call_filter_cn_neutral_sex_chrom(self):
        """--filter cn respects expected CN on sex chromosomes (male sample)."""
        segarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1", "chr1", "chrX", "chrX"],
                    "start": [0, 1000, 0, 1000],
                    "end": [1000, 2000, 1000, 2000],
                    "gene": ["A", "B", "C", "D"],
                    "log2": [0.0, 0.01, 0.0, 0.05],
                    "weight": [1.0, 1.0, 1.0, 1.0],
                    "probes": [10, 10, 10, 10],
                }
            )
        )
        # Male sample: expected CN=2 on autosomes, CN=1 on chrX
        result = commands.do_call(
            segarr,
            method="threshold",
            is_haploid_x_reference=True,
            is_sample_female=False,
            filters=["cn"],
        )
        # chr1 A,B are cn=2 (neutral for autosome) — merged
        # chrX C,D are cn=1 (neutral for male X) — merged
        self.assertEqual(len(result), 2)

        # Female sample: expected CN=2 on chrX too
        result_f = commands.do_call(
            segarr,
            method="threshold",
            is_haploid_x_reference=True,
            is_sample_female=True,
            filters=["cn"],
        )
        # chr1 A,B are cn=2 (neutral) — merged
        # chrX C,D are cn=1 (NOT neutral for female, expected 2) — kept separate
        self.assertEqual(len(result_f), 3)

    def test_call_log2_ratios(self):
        cnarr = cnvlib.read("formats/par-reference.grch38.cnn")
        ploidy = 2
        purity = 0.8
        is_haploid_x_reference = True
        is_sample_female = False
        diploid_parx_genome = None
        absolutes = call.absolute_clonal(
            cnarr,
            ploidy,
            purity,
            is_haploid_x_reference,
            diploid_parx_genome,
            is_sample_female,
        )
        ratios = call.log2_ratios(
            cnarr, absolutes, ploidy, is_haploid_x_reference, diploid_parx_genome
        )
        ratios_dprx = call.log2_ratios(
            cnarr, absolutes, ploidy, is_haploid_x_reference, "grch38"
        )
        self.assertEqual(len(ratios), len(cnarr))
        self.assertEqual(len(ratios_dprx), len(cnarr))
        self.assertEqual(ratios[26], ratios_dprx[26])
        self.assertEqual(ratios[26], ratios_dprx[27] + 1)

    def test_call_sex(self):
        """Test each 'call' method on allosomes."""
        for (
            fname,
            sample_is_f,
            ref_is_m,
            chr1_expect,
            chrx_expect,
            chry_expect,
            chr1_cn,
            chrx_cn,
            chry_cn,
        ) in (
            ("formats/f-on-f.cns", True, False, 0, 0, None, 2, 2, None),
            ("formats/f-on-m.cns", True, True, 0.585, 1, None, 3, 2, None),
            ("formats/m-on-f.cns", False, False, 0, -1, 0, 2, 1, 1),
            ("formats/m-on-m.cns", False, True, 0, 0, 0, 2, 1, 1),
        ):
            cns = cnvlib.read(fname)
            chr1_idx = cns.chromosome == "chr1"
            chrx_idx = cns.chromosome == "chrX"
            chry_idx = cns.chromosome == "chrY"

            def test_chrom_means(segments):
                self.assertEqual(chr1_cn, segments["cn"][chr1_idx].mean())
                self.assertAlmostEqual(
                    chr1_expect, segments["log2"][chr1_idx].mean(), 0
                )
                self.assertEqual(chrx_cn, segments["cn"][chrx_idx].mean())
                self.assertAlmostEqual(
                    chrx_expect, segments["log2"][chrx_idx].mean(), 0
                )
                if not sample_is_f:
                    self.assertEqual(chry_cn, segments["cn"][chry_idx].mean())
                    self.assertAlmostEqual(
                        chry_expect, segments["log2"][chry_idx].mean(), 0
                    )

            # Call threshold
            cns_thresh = commands.do_call(
                cns,
                None,
                "threshold",
                is_haploid_x_reference=ref_is_m,
                is_sample_female=sample_is_f,
            )
            test_chrom_means(cns_thresh)
            # Call clonal pure
            cns_clone = commands.do_call(
                cns,
                None,
                "clonal",
                is_haploid_x_reference=ref_is_m,
                is_sample_female=sample_is_f,
            )
            test_chrom_means(cns_clone)
            # Call clonal barely-mixed
            cns_p99 = commands.do_call(
                cns,
                None,
                "clonal",
                purity=0.99,
                is_haploid_x_reference=ref_is_m,
                is_sample_female=sample_is_f,
            )
            test_chrom_means(cns_p99)

    def test_call_various_abs_ref_exp_methods(self):
        cnarr = cnvlib.read("formats/par-reference.grch38.cnn")

        def _run(is_haploid_x_reference, is_sample_female, diploid_parx_genome=None):
            ploidy = 2
            purity = 0.8
            abs_df = call.absolute_dataframe(
                cnarr,
                ploidy,
                purity,
                is_haploid_x_reference,
                diploid_parx_genome,
                is_sample_female,
            )
            abs_ref = call.absolute_reference(
                cnarr, ploidy, diploid_parx_genome, is_haploid_x_reference
            )
            abs_exp = call.absolute_expect(
                cnarr, ploidy, diploid_parx_genome, is_sample_female
            )
            abs_clonal = call.absolute_clonal(
                cnarr,
                ploidy,
                purity,
                is_haploid_x_reference,
                diploid_parx_genome,
                is_sample_female,
            )
            return abs_df, abs_ref, abs_exp, abs_clonal

        def _assert_abs_df(iloc, abs_df, ref_copies, exp_copies):
            self.assertTrue("reference" in abs_df.columns)
            self.assertTrue("expect" in abs_df.columns)
            r = abs_df.iloc[iloc]
            self.assertEqual(r.reference, ref_copies)
            self.assertEqual(r.expect, exp_copies)

        def _assert_abs_copies(i, abs_values, copies):
            self.assertEqual(abs_values[i], copies)

        def _assert_abs_clonal(i, abs_clonal, value):
            self.assertAlmostEqual(abs_clonal[i], value, 5)

        def _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal):
            i = 0
            _assert_abs_df(i, abs_df, 2, 2)
            _assert_abs_copies(i, abs_ref, 2)
            _assert_abs_copies(i, abs_exp, 2)
            _assert_abs_clonal(i, abs_clonal, 0.26708)

        def _assert_chrx_par(
            abs_df, abs_ref, abs_exp, abs_clonal, ref_copies, exp_copies, clonal_copies
        ):
            i = 13
            _assert_abs_df(i, abs_df, ref_copies, exp_copies)
            _assert_abs_copies(i, abs_ref, ref_copies)
            _assert_abs_copies(i, abs_exp, exp_copies)
            _assert_abs_clonal(i, abs_clonal, clonal_copies)

        def _assert_chrx_non_par(
            abs_df, abs_ref, abs_clonal, abs_exp, ref_copies, exp_copies, clonal_copies
        ):
            i = 21
            _assert_abs_df(i, abs_df, ref_copies, exp_copies)
            _assert_abs_copies(i, abs_ref, ref_copies)
            _assert_abs_copies(i, abs_exp, exp_copies)
            _assert_abs_clonal(i, abs_clonal, clonal_copies)

        def _assert_chry_par(
            abs_df, abs_ref, abs_exp, abs_clonal, ref_copies, exp_copies, clonal_copies
        ):
            i = 36
            _assert_abs_df(i, abs_df, ref_copies, exp_copies)
            _assert_abs_copies(i, abs_ref, ref_copies)
            _assert_abs_copies(i, abs_exp, exp_copies)
            _assert_abs_clonal(i, abs_clonal, clonal_copies)

        def _assert_chry_non_par(
            abs_df, abs_ref, abs_exp, abs_clonal, ref_copies, exp_copies, clonal_copies
        ):
            i = 40
            _assert_abs_df(i, abs_df, ref_copies, exp_copies)
            _assert_abs_copies(i, abs_ref, ref_copies)
            _assert_abs_copies(i, abs_exp, exp_copies)
            _assert_abs_clonal(i, abs_clonal, clonal_copies)

        is_haploid_x_reference = True
        is_female_sample = True
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 2, 1.59225)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 1, 2, 1.34001)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 0, 2.04202)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 0, 0.70083)
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample, "grch38"
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 2, 2, 3.68449)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 1, 2, 1.34001)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 0, 0, 0.0)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 0, 0.70083)

        is_haploid_x_reference = True
        is_female_sample = False
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 1, 1.84225)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 1, 1, 1.59001)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 1, 1.79202)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 1, 0.45083)
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample, "grch38"
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 2, 2, 3.68449)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 1, 1, 1.59001)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 0, 0, 0.0)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 1, 0.45083)

        is_haploid_x_reference = False
        is_female_sample = True
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 2, 2, 3.68449)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 2, 2, 3.18002)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 0, 2.04202)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 0, 0.70083)
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample, "grch38"
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 2, 2, 3.684493)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 2, 2, 3.18002)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 0, 0, 0.0)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 0, 0.70083)

        is_haploid_x_reference = False
        is_female_sample = False
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 2, 1, 3.93449)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 2, 1, 3.43002)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 1, 1.79202)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 1, 0.45083)
        abs_df, abs_ref, abs_exp, abs_clonal = _run(
            is_haploid_x_reference, is_female_sample, "grch38"
        )
        _assert_chr1(abs_df, abs_ref, abs_exp, abs_clonal)
        _assert_chrx_par(abs_df, abs_ref, abs_exp, abs_clonal, 2, 2, 3.68449)
        _assert_chrx_non_par(abs_df, abs_ref, abs_clonal, abs_exp, 2, 1, 3.43002)
        _assert_chry_par(abs_df, abs_ref, abs_exp, abs_clonal, 0, 0, 0.0)
        _assert_chry_non_par(abs_df, abs_ref, abs_exp, abs_clonal, 1, 1, 0.45083)

    @staticmethod
    def _deep_deletion_cns():
        """Segments including deep deletions: a true homozygous deletion
        (log2=-3) and a zero-coverage sentinel segment (log2=NULL_LOG2_COVERAGE).
        """
        log2 = [0.0, -3.0, params.NULL_LOG2_COVERAGE, 0.585]
        n = len(log2)
        return cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(0, n * 1000, 1000),
                    "end": np.arange(1000, n * 1000 + 1000, 1000),
                    "gene": ["-"] * n,
                    "log2": log2,
                    "probes": [10] * n,
                    "weight": [1.0] * n,
                }
            )
        )

    def test_call_clonal_no_negative_cn(self):
        """Clonal calling with impure samples never emits negative copy number.

        Regression for #503/#516: with purity < 1, the purity-rescale formula
        n = (r*2^log2 - x*(1-p)) / p extrapolates to negative absolute copies
        for deeply deleted (or zero-coverage sentinel) segments. Absolute copy
        number is physically >= 0, so it must be floored at 0.
        """
        cns = self._deep_deletion_cns()
        for tumor_purity in (0.3, 0.5, 0.65, 0.9):
            with self.subTest(purity=tumor_purity):
                called = commands.do_call(
                    cns,
                    None,
                    "clonal",
                    purity=tumor_purity,
                    is_haploid_x_reference=True,
                    is_sample_female=True,
                )
                self.assertTrue(
                    (called["cn"] >= 0).all(),
                    f"negative cn at purity={tumor_purity}: {called['cn'].tolist()}",
                )
                # Deep-deletion and sentinel segments should call CN 0
                self.assertEqual(called["cn"].iloc[1], 0)
                self.assertEqual(called["cn"].iloc[2], 0)

    def test_log2_ratio_to_absolute_floored_at_zero(self):
        """_log2_ratio_to_absolute never returns a negative absolute (#503)."""
        # Impure path: deep deletion below the (1-p) contamination floor
        self.assertEqual(call._log2_ratio_to_absolute(-20.0, 2, 2, purity=0.5), 0.0)
        self.assertEqual(call._log2_ratio_to_absolute(-3.0, 2, 2, purity=0.5), 0.0)
        # A real single-copy gain is unaffected
        self.assertAlmostEqual(
            call._log2_ratio_to_absolute(0.585, 2, 2, purity=1.0), 3.0, places=2
        )
        # A NaN log2 (malformed input) propagates rather than being silently
        # floored to 0 -- a false homozygous-deletion call would be worse.
        self.assertTrue(
            np.isnan(call._log2_ratio_to_absolute(np.nan, 2, 2, purity=0.5))
        )


class LoadHetSnpsTests(unittest.TestCase):
    """Tests for cmdutil.load_het_snps and the _warn_if_baf_input_suspicious helper."""

    def test_helper_warns_on_empty(self):
        """Empty heterozygous result emits a clear warning."""
        with self.assertLogs(level="WARNING") as ctx:
            cmdutil._warn_if_baf_input_suspicious(None)
        self.assertTrue(
            any("No heterozygous variants" in msg for msg in ctx.output),
            ctx.output,
        )

    def test_helper_warns_on_skewed_distribution(self):
        """Median alt_freq far from 0.5 emits a warning."""
        # 100 variants all near 1.0 (e.g. mistakenly homozygous-alt or somatic)
        alt_freqs = pd.Series([0.95] * 100)
        with self.assertLogs(level="WARNING") as ctx:
            cmdutil._warn_if_baf_input_suspicious(alt_freqs)
        self.assertTrue(
            any("Median allele frequency" in msg for msg in ctx.output),
            ctx.output,
        )

    def test_helper_silent_on_balanced_distribution(self):
        """Median alt_freq near 0.5 emits no warning."""
        rng = np.random.default_rng(0)
        # 100 het SNPs with alt_freq centered on 0.5 (binomial-ish noise)
        alt_freqs = pd.Series(rng.normal(0.5, 0.05, 100).clip(0, 1))
        with self.assertNoLogs(level="WARNING"):
            cmdutil._warn_if_baf_input_suspicious(alt_freqs)

    def test_helper_silent_below_min_count(self):
        """Distribution check is skipped for small variant sets to avoid noise."""
        # Skewed but only 10 variants -- no warning (could be a small panel)
        alt_freqs = pd.Series([0.95] * 10)
        with self.assertNoLogs(level="WARNING"):
            cmdutil._warn_if_baf_input_suspicious(alt_freqs)

    def test_load_het_snps_warns_on_empty_vcf(self):
        """Loading a VCF with no records triggers the empty-result warning."""
        with self.assertLogs(level="WARNING") as ctx:
            varr = commands.load_het_snps("formats/blank.vcf", None, None, 1, None)
        self.assertEqual(len(varr), 0)
        self.assertTrue(
            any("No heterozygous variants" in msg for msg in ctx.output),
            ctx.output,
        )

    def test_load_het_snps_silent_on_real_vcf(self):
        """A legitimate test VCF with paired sample IDs doesn't trigger warnings."""
        with self.assertNoLogs(level="WARNING"):
            varr = commands.load_het_snps(
                "formats/na12878_na12882_mix.vcf", "NA12878", "NA12882", 15, None
            )
        self.assertGreater(len(varr), 50)

    def test_load_het_snps_warns_on_missing_sample_ids(self):
        """Without -i/-n, the wrong sample is used; the distribution warning fires."""
        with self.assertLogs(level="WARNING") as ctx:
            commands.load_het_snps(
                "formats/na12878_na12882_mix.vcf", None, None, 15, None
            )
        self.assertTrue(
            any("Median allele frequency" in msg for msg in ctx.output),
            ctx.output,
        )
