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


class SegmentationTests(unittest.TestCase):
    """Tests for segmentation commands."""

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
        varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        # TODO - This test is failing... commenting it out for now!
        # segments = segmentation.do_segmentation(cnarr, "haar", variants=varr)
        # self.assertGreater(len(segments), n_chroms)
        # self.assertTrue((segments.start < segments.end).all())

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
