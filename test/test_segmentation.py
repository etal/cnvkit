#!/usr/bin/env python
"""Tests for segmentation commands (Python-side integration)."""

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


class SegmentationTests(unittest.TestCase):
    """Tests for segmentation commands."""

    def test_segment_warns_on_missing_sample_id(self):
        """A CopyNumArray with no sample_id emits one warning and falls back to
        the literal 'None' as the segment ID. The CLI never reaches this --
        read() derives sample_id from the input filename -- but the raw
        in-memory API can construct an array without one.
        """
        cnarr = cnvlib.read("formats/amplicon.cnr")
        cnarr.meta.pop("sample_id", None)
        self.assertIsNone(cnarr.sample_id)
        with self.assertLogs(level="WARNING") as cm:
            segmentation.do_segmentation(cnarr, "haar")
        self.assertTrue(any("sample_id" in line for line in cm.output))

    def test_segment_no_warning_when_sample_id_present(self):
        """The missing-sample_id warning stays silent when sample_id is set."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        self.assertIsNotNone(cnarr.sample_id)
        with self.assertLogs(level="INFO") as cm:
            segmentation.do_segmentation(cnarr, "haar")
        self.assertFalse(any("no sample_id" in line.lower() for line in cm.output))

    def test_segment(self):
        """The 'segment' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        n_chroms = cnarr.chromosome.nunique()
        # NB: R methods are in another script; haar is pure-Python
        segments = segmentation.do_segmentation(cnarr, "haar")
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())
        segments = segmentation.do_segmentation(
            cnarr, "haar", threshold=0.0001, skip_low=True
        )
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())
        # haar segmentation with variants unions depth+BAF breakpoints
        # (cnvkit-ugh); see test_haar_vcf_detects_copy_neutral_loh for the LOH
        # behavior. Here just confirm it runs and yields valid segments.
        varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        segments = segmentation.do_segmentation(cnarr, "haar", variants=varr)
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())

    @pytest.mark.slow
    def test_segment_hmm(self):
        """The 'segment' command with HMM methods."""
        # Test all HMM method variants on one file
        cnarr = cnvlib.read("formats/amplicon.cnr")
        n_chroms = cnarr.chromosome.nunique()
        # NB: R methods are in another script; haar is pure-Python
        segments = segmentation.do_segmentation(cnarr, "hmm")
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())
        segments = segmentation.do_segmentation(cnarr, "hmm-tumor", skip_low=True)
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())
        segments = segmentation.do_segmentation(cnarr, "hmm-germline")
        self.assertGreater(len(segments), n_chroms)
        self.assertTrue((segments.start < segments.end).all())
        varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        segments = segmentation.do_segmentation(cnarr, "hmm", variants=varr)
        self.assertGreater(len(segments), n_chroms)
        # Verify default HMM also works on a different dataset
        cnarr2 = cnvlib.read("formats/p2-20_1.cnr")
        n_chroms2 = cnarr2.chromosome.nunique()
        segments = segmentation.do_segmentation(cnarr2, "hmm")
        self.assertGreater(len(segments), n_chroms2)
        self.assertTrue((segments.start < segments.end).all())

    @pytest.mark.slow
    def test_segment_parallel(self):
        """The 'segment' command, in parallel."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        psegments = segmentation.do_segmentation(cnarr, "haar", processes=2)
        ssegments = segmentation.do_segmentation(cnarr, "haar", processes=1)
        self.assertEqual(psegments.data.shape, ssegments.data.shape)
        self.assertEqual(len(psegments.meta), len(ssegments.meta))
        # Parallel and serial must agree on segment boundaries and values, not
        # just shape -- a chromosome-ordering or result-assembly bug would slip
        # past a shape-only check.
        psorted = psegments.data.sort_values(["chromosome", "start"]).reset_index(
            drop=True
        )
        ssorted = ssegments.data.sort_values(["chromosome", "start"]).reset_index(
            drop=True
        )
        self.assertEqual(list(psorted["chromosome"]), list(ssorted["chromosome"]))
        np.testing.assert_array_equal(
            psorted["start"].to_numpy(), ssorted["start"].to_numpy()
        )
        np.testing.assert_array_equal(
            psorted["end"].to_numpy(), ssorted["end"].to_numpy()
        )
        np.testing.assert_allclose(
            psorted["log2"].to_numpy(), ssorted["log2"].to_numpy(), atol=1e-9
        )

    def test_segment_empty_input(self):
        """Test segmentation with empty CNR input (issue #970)."""
        # Create an empty CNA with proper structure (header only)
        empty_data = pd.DataFrame(
            columns=["chromosome", "start", "end", "gene", "log2"]
        )
        empty_cnarr = cnvlib.cnary.CopyNumArray(empty_data, {"sample_id": "test"})

        # Test with serial processing
        segments = segmentation.do_segmentation(empty_cnarr, "haar", processes=1)
        self.assertEqual(len(segments), 0)
        self.assertListEqual(
            list(segments.data.columns), list(empty_cnarr.data.columns)
        )

        # Test with parallel processing
        psegments = segmentation.do_segmentation(empty_cnarr, "haar", processes=2)
        self.assertEqual(len(psegments), 0)

        # Test with save_dataframe=True
        segments_df, rstr = segmentation.do_segmentation(
            empty_cnarr, "haar", processes=1, save_dataframe=True
        )
        self.assertEqual(len(segments_df), 0)
        self.assertEqual(rstr, "")

    def test_threshold_zero_not_overridden(self):
        """threshold=0.0 is honored, not treated as 'unset' (scq.4).

        'if not threshold' replaced the valid value 0.0 with the method default;
        'if threshold is None' keeps it.
        """
        cnarr = cnvlib.read("formats/amplicon.cnr")
        with self.assertLogs(level="INFO") as cm:
            segmentation.do_segmentation(cnarr, "haar", threshold=0.0, processes=1)
        msg = " ".join(cm.output)
        self.assertIn("significance threshold 0.0,", msg)
        self.assertNotIn("0.0001", msg)

    def test_haar_vcf_detects_copy_neutral_loh(self):
        """`-m haar --vcf` unions depth+BAF breakpoints -> catches copy-neutral
        LOH (flat depth, BAF shift), and adds a 'baf' column (cnvkit-ugh)."""
        n = 120
        rng = np.random.default_rng(0)
        cnarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(n) * 1000,
                    "end": np.arange(n) * 1000 + 1000,
                    "gene": "-",
                    "log2": rng.normal(0.0, 0.05, n),  # flat depth, no CN change
                    "weight": np.ones(n),
                }
            )
        )
        # Het SNPs: balanced (minor/depth=0.5) first half, LOH (0.3) second half
        alt = np.concatenate([np.full(60, 15.0), np.full(60, 9.0)])
        varr = vary.VariantArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(n) * 1000 + 500,
                    "end": np.arange(n) * 1000 + 501,
                    "ref": ["A"] * n,
                    "alt": ["G"] * n,
                    "zygosity": np.full(n, 0.5),
                    "alt_count": alt,
                    "depth": np.full(n, 30.0),
                    "alt_freq": alt / 30.0,
                }
            )
        )
        # Depth-only: flat -> no BAF breakpoint, no 'baf' column
        depth_only = segmentation.do_segmentation(cnarr, "haar")
        self.assertNotIn("baf", depth_only.data.columns)
        # Joint: the BAF shift at bin 60 must introduce a breakpoint there
        joint = segmentation.do_segmentation(cnarr, "haar", variants=varr)
        self.assertIn("baf", joint.data.columns)
        self.assertGreater(len(joint), len(depth_only))
        boundaries = set(joint.start.tolist()) | set(joint.end.tolist())
        self.assertTrue(
            any(55000 <= b <= 65000 for b in boundaries),
            f"expected a breakpoint near the BAF shift (~60000); got {sorted(boundaries)}",
        )

    def test_haar_vcf_without_allele_info_falls_back(self):
        """`-m haar --vcf` on a VCF lacking AF and AD/DP must not crash; it
        falls back to depth-only segmentation (cnvkit-ugh)."""
        n = 30
        cnarr = cnary.CopyNumArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * n,
                    "start": np.arange(n) * 1000,
                    "end": np.arange(n) * 1000 + 1000,
                    "gene": "-",
                    "log2": np.zeros(n),
                    "weight": np.ones(n),
                }
            )
        )
        # Genotype-only VCF: het zygosity, but no alt_freq / alt_count / depth
        varr = vary.VariantArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 10,
                    "start": np.arange(10) * 100,
                    "end": np.arange(10) * 100 + 1,
                    "ref": ["A"] * 10,
                    "alt": ["G"] * 10,
                    "zygosity": np.full(10, 0.5),
                }
            )
        )
        seg = segmentation.do_segmentation(cnarr, "haar", variants=varr)
        self.assertGreater(len(seg), 0)
        self.assertTrue((seg.start < seg.end).all())


class TransferFieldsTests(unittest.TestCase):
    """Tests for transfer_fields and do_segmentation NaN handling."""

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
                    "start": [
                        54246732,
                        54800000,
                        55018770,
                        55019096,
                        55032092,
                        55088166,
                    ],
                    "end": [54300000, 54900000, 55019096, 55019423, 55032193, 55088469],
                    "gene": ["VSTM2A", "SEC61G", "EGFR", "EGFR", "EGFR", "EGFR"],
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

    def test_do_segmentation_drops_inf_log2(self):
        """do_segmentation tolerates ±inf-log2 bins (gh#508, sibling to #881).

        The #881 fix dropped NaN-log2 bins before drop_outliers' Savitzky-Golay
        smoother because scipy's lstsq rejects non-finite input ("array must
        not contain infs or NaNs"). But pandas ``.isna()`` does NOT catch
        ±inf, so degenerate flat-reference WGS data (gh#508) -- where some
        bins can land at ±inf after reference subtraction -- still crashed
        on the same path. Broadening the prefilter from ``.isna()`` to
        ``~np.isfinite`` covers both.
        """
        n = 120  # > drop_outliers' width (50) so savgol actually runs
        rng = np.random.default_rng(0)
        log2 = rng.normal(0, 0.2, n)
        log2[5] = np.inf  # inside the leading savgol edge window
        log2[60] = -np.inf
        log2[80] = np.nan  # confirm NaN is still handled alongside ±inf
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
        # Must not raise, and the segment log2 must be finite (no inf/nan leaks).
        cns = segmentation.do_segmentation(cnarr, "haar", processes=1)
        self.assertGreater(len(cns), 0)
        self.assertTrue(np.isfinite(cns["log2"].to_numpy()).all())
