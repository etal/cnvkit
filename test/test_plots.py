#!/usr/bin/env python
"""Tests for plotting commands: scatter, heatmap, diagram."""

import ast
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


class PlotTests(unittest.TestCase):
    """Smoke tests for plotting commands."""

    def test_scatter(self):
        """The 'scatter' command."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        fig = scatter.do_scatter(cnarr, segarr)
        self.assertIsNotNone(fig)
        # With a gene zoom
        fig = scatter.do_scatter(cnarr, segarr, show_gene="BRAF")
        self.assertIsNotNone(fig)

    def test_scatter_genome_y_floor(self):
        """Genome-wide autoscale floors y_min so a single deep homozygous
        deletion can't compress the whole plot (#385)."""
        # lazy: defer matplotlib import to keep headless test collection fast
        from matplotlib import pyplot  # noqa: PLC0415

        probes = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 200, "A", 0.0],
                ["chr1", 200, 300, "A", -0.1],
                ["chr1", 300, 400, "A", -12.0],
                ["chr2", 100, 200, "B", 0.05],
                ["chr2", 200, 300, "B", -0.2],
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        segments = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 400, "A", -12.0],
                ["chr2", 100, 300, "B", 0.0],
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        segments.data["probes"] = [3, 2]
        _fig, ax = pyplot.subplots()
        scatter.cnv_on_genome(ax, probes, segments)
        y_lo, _y_hi = ax.get_ylim()
        pyplot.close(_fig)
        self.assertGreaterEqual(
            y_lo,
            scatter.AUTO_Y_MIN_FLOOR,
            "Deep deletion must not pull the genome-wide y-axis below the floor",
        )
        # An explicit --y-min still overrides the floor
        _fig, ax = pyplot.subplots()
        scatter.cnv_on_genome(ax, probes, segments, y_min=-15.0, y_max=2.0)
        y_lo, _y_hi = ax.get_ylim()
        pyplot.close(_fig)
        self.assertAlmostEqual(y_lo, -15.0)

    def test_scatter_show_snvs_extras(self):
        """`scatter` overlays LOH and somatic markers on the VAF panel and
        leaves the segment-VAF trend driven solely by the het subset (#290).

        The structural-correctness invariant: introducing the overlays must
        NOT pull the segment-mean VAF toward the homozygous-allele extremes,
        because clinical interpretation of segment VAFs depends on the trend
        being the mean of the het BAFs only. We verify by monkey-patching
        ``get_segment_vafs`` to capture its first positional argument (the
        VariantArray driving the trend) and asserting it's the het subset
        with no LOH rows mixed in -- the test catches a regression where
        ``snv_on_chromosome`` accidentally concatenates LOH into the trend.
        """
        # lazy: defer matplotlib import to keep headless test collection fast
        from matplotlib import pyplot  # noqa: PLC0415

        het = vary.VariantArray.from_rows(
            [
                ("chr1", 100, 101, "A", "G", False, 0.5, 0.5, 100, 50),
                ("chr1", 200, 201, "A", "G", False, 0.5, 0.55, 100, 55),
                ("chr1", 300, 301, "A", "G", False, 0.5, 0.45, 100, 45),
            ],
            columns=[
                "chromosome",
                "start",
                "end",
                "ref",
                "alt",
                "somatic",
                "zygosity",
                "alt_freq",
                "depth",
                "alt_count",
            ],
        )
        loh = vary.VariantArray.from_rows(
            [
                ("chr1", 150, 151, "A", "G", False, 1.0, 0.95, 100, 95),
                ("chr1", 250, 251, "A", "G", False, 0.0, 0.05, 100, 5),
            ],
            columns=[
                "chromosome",
                "start",
                "end",
                "ref",
                "alt",
                "somatic",
                "zygosity",
                "alt_freq",
                "depth",
                "alt_count",
            ],
        )
        seg = cnary.CopyNumArray.from_rows(
            [("chr1", 50, 400, "A", -0.5)],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        # Monkey-patch get_segment_vafs to capture the variants array it sees.
        captured = []
        orig_get_segment_vafs = scatter.get_segment_vafs

        def recording_get_segment_vafs(variants, segments):
            captured.append(variants)
            return orig_get_segment_vafs(variants, segments)

        scatter.get_segment_vafs = recording_get_segment_vafs
        try:
            _fig, ax = pyplot.subplots()
            scatter.snv_on_chromosome(
                ax,
                het,
                seg,
                [],
                do_trend=True,
                by_bin=False,
                segment_color="darkorange",
                loh_variants=loh,
            )
            pyplot.close(_fig)
        finally:
            scatter.get_segment_vafs = orig_get_segment_vafs

        # snv_on_chromosome must drive the trend from the het subset only.
        self.assertEqual(len(captured), 1)
        trend_input = captured[0]
        self.assertEqual(len(trend_input), len(het))
        # The LOH rows' VAF extremes (0.95, 0.05) must NOT be among the
        # values the trend was computed from. If they were, the test caught
        # a regression where the overlay leaked into the BAF math.
        for loh_freq in loh["alt_freq"].tolist():
            self.assertNotIn(loh_freq, trend_input["alt_freq"].tolist())
        # And as a sanity check on the test's premise: a "broken" trend
        # input that concatenates LOH would produce a clearly different
        # segment mean.
        broken_input = het.concat([loh])
        broken_trend = [
            v_freq for _seg, v_freq in orig_get_segment_vafs(broken_input, seg)
        ]
        correct_trend = [
            v_freq for _seg, v_freq in orig_get_segment_vafs(trend_input, seg)
        ]
        self.assertNotEqual(broken_trend, correct_trend)

    def test_heatmap(self):
        """The 'heatmap' command."""
        cnarrs = [cnvlib.read("formats/amplicon.cnr")]
        ax = heatmap.do_heatmap(cnarrs)
        self.assertIsNotNone(ax)
        # With desaturation
        ax = heatmap.do_heatmap(cnarrs, do_desaturate=True)
        self.assertIsNotNone(ax)


class GeneCoordsTests(unittest.TestCase):
    """Tests for plots.gene_coords_by_name gene-label selection."""

    def test_gene_coords_by_name(self):
        """`-g` labels only requested genes, not co-binned neighbors (#458).

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
        for _s, _e, label in plots.gene_coords_by_name(cnarr, ["ERBB2", "MIR4728"])[
            "chr17"
        ]:
            names_seen.update(label.split(","))
        self.assertEqual(names_seen, {"ERBB2", "MIR4728"})


class DiagramGeneLabelTests(unittest.TestCase):
    """`diagram` --gene and directional --threshold-low/-high (#248)."""

    @staticmethod
    def _segments():
        """A .cns with one clear gain, one clear loss, and a co-binned label."""
        seg = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 200, "GAINER", 1.0],
                ["chr2", 100, 200, "LOSER", -1.0],
                ["chr3", 100, 200, "FLAT", 0.1],
                ["chr17", 100, 200, "NF1,ERBB2,ERBB2,MIR4728", 0.8],
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        seg.data["probes"] = [5, 5, 5, 5]
        return seg

    def test_symmetric_threshold_unchanged(self):
        """Default symmetric threshold labels both gains and losses."""
        seg = self._segments()
        labels = diagram._get_gene_labels(seg, None, True, -0.5, 0.5, 3)
        self.assertIn("GAINER", labels)
        self.assertIn("LOSER", labels)
        self.assertNotIn("FLAT", labels)

    def test_threshold_high_labels_only_gains(self):
        """--threshold-high (loss side disabled) labels only gains."""
        seg = self._segments()
        labels = diagram._get_gene_labels(seg, None, True, None, 0.5, 3)
        self.assertIn("GAINER", labels)
        self.assertNotIn("LOSER", labels)

    def test_threshold_low_labels_only_losses(self):
        """--threshold-low (gain side disabled) labels only losses."""
        seg = self._segments()
        labels = diagram._get_gene_labels(seg, None, True, -0.5, None, 3)
        self.assertIn("LOSER", labels)
        self.assertNotIn("GAINER", labels)

    def test_threshold_low_suppresses_deletions(self):
        """A very low loss cutoff (e.g. -25) suppresses all loss labels while a
        gain cutoff still surfaces gains (maintainer's 2018 example)."""
        seg = self._segments()
        labels = diagram._get_gene_labels(seg, None, True, -25.0, 0.5, 3)
        self.assertIn("GAINER", labels)
        self.assertNotIn("LOSER", labels)

    def test_min_probes_filters(self):
        """Below-min-probes segments are not labeled regardless of threshold."""
        seg = self._segments()
        seg.data["probes"] = [1, 1, 1, 1]
        labels = diagram._get_gene_labels(seg, None, True, -0.5, 0.5, 3)
        self.assertEqual(labels, [])

    def test_gene_feature_label_default(self):
        """Without -g, the comma-joined gene string is expanded verbatim."""
        self.assertEqual(
            diagram._gene_feature_label("NF1,ERBB2,ERBB2,MIR4728", None),
            "NF1, ERBB2, ERBB2, MIR4728",
        )
        self.assertEqual(diagram._gene_feature_label("GAINER", None), "GAINER")

    def test_gene_feature_label_restricts_to_requested(self):
        """-g shows only requested genes, dropping co-binned neighbors and dups."""
        # ERBB2 requested: NF1/MIR4728 neighbors dropped; duplicate ERBB2 collapsed
        self.assertEqual(
            diagram._gene_feature_label("NF1,ERBB2,ERBB2,MIR4728", {"ERBB2"}),
            "ERBB2",
        )
        # Multiple requested genes sharing the bin all appear, in order
        self.assertEqual(
            diagram._gene_feature_label("NF1,ERBB2,MIR4728", {"ERBB2", "MIR4728"}),
            "ERBB2, MIR4728",
        )
        # A segment with no requested gene yields no label
        self.assertIsNone(diagram._gene_feature_label("STARD3", {"ERBB2"}))
        # Whitespace around comma-separated names is tolerated when matching
        self.assertEqual(
            diagram._gene_feature_label("NF1, ERBB2", {"ERBB2"}),
            "ERBB2",
        )

    def test_create_diagram_smoke_with_flags(self):
        """End-to-end: the new flags drive create_diagram and write a PDF."""
        seg = self._segments()
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "diagram.pdf")
            result = diagram.create_diagram(
                None,
                seg,
                0.5,
                3,
                out,
                title="Test sample",
                threshold_high=0.5,
                gene_names=["ERBB2"],
            )
            self.assertTrue(os.path.exists(result))

    def test_cmd_diagram_threshold_mutually_exclusive(self):
        """Explicit -t with a directional threshold is rejected up front."""
        args = commands.parse_args(
            ["diagram", "formats/amplicon.cns", "-t", "0.5", "--threshold-low", "-0.5"]
        )
        with self.assertRaises(ValueError):
            args.func(args)

    def test_cmd_diagram_plumbs_new_args_to_create_diagram(self):
        """`_cmd_diagram` forwards gene_names and directional thresholds.

        AST-walk plumbing check (no fixtures): assert the create_diagram call
        site passes the new keyword arguments, so a dropped kwarg fails fast.
        """
        tree = ast.parse(inspect.getsource(commands._cmd_diagram))
        kwargs_seen = set()
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "create_diagram"
            ):
                kwargs_seen = {kw.arg for kw in node.keywords}
        self.assertLessEqual(
            {"threshold_low", "threshold_high", "gene_names"}, kwargs_seen
        )


class DiagramCoordinateTests(unittest.TestCase):
    """diagram renders reverse-oriented intervals (start > end) gracefully.

    Reverse-direction PCR primers are sometimes saved in the input BED with
    start > end; CNVkit should treat the interval as spanning [min, max] rather
    than crashing Biopython's renderer, which asserts start <= end <= length.
    """

    def test_feature_span_normalizes_orientation(self):
        # Forward interval: 1-based start -> 0-based half-open, unchanged
        self.assertEqual(diagram._feature_span(100, 200, 300), (99, 200))
        # Reversed interval (start > end): normalized to [min-1, max)
        self.assertEqual(diagram._feature_span(400, 300, 400), (299, 400))
        # Out of range (beyond chromosome) is rejected
        self.assertIsNone(diagram._feature_span(100, 500, 300))
        # start == 0 -> lo == -1 is rejected, preserving prior lower-bound guard
        self.assertIsNone(diagram._feature_span(0, 50, 300))

    def test_chromosome_sizes_counts_reversed_intervals(self):
        seg = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 200, "FWD", 0.0],
                ["chr1", 400, 300, "REV", 0.0],  # reversed: rightmost pos is 400
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        sizes = plots.chromosome_sizes(seg)
        self.assertEqual(sizes["chr1"], 400)

    def test_create_diagram_renders_reversed_interval(self):
        seg = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 200, "FWD", 1.0],
                ["chr1", 400, 300, "REV", 1.0],  # reverse-primer style
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        seg.data["probes"] = [5, 5]
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "diagram.pdf")
            # Must not raise the Biopython start<=end<=length assertion
            result = diagram.create_diagram(None, seg, 0.5, 3, out, title="t")
            self.assertTrue(os.path.exists(result))

    def test_create_diagram_amplicon_fixture(self):
        """Regression: the amplicon fixture has several reverse-oriented
        segments (e.g. chr2 212578209-212576985) that previously crashed the
        renderer. It must now produce a diagram."""
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "diagram.pdf")
            result = diagram.create_diagram(cnarr, segarr, 0.5, 3, out)
            self.assertTrue(os.path.exists(result))
