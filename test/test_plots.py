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

    Such a row cannot reach here from a file -- `skgenome.tabio.read`
    reverses it out of a BED and refuses it in CNVkit's own formats -- so
    these build their rows in memory, which is the path that remains.
    Biopython's renderer asserts start <= end <= length, and used to crash.
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


class ByBinCoordinateTests(unittest.TestCase):
    """Both --by-bin lookups refuse rows their axis cannot represent.

    The rationale lives with the code, in `plots.update_binwise_positions`.
    """

    @staticmethod
    def _bins(rows):
        return cnary.CopyNumArray.from_rows(
            rows, columns=["chromosome", "start", "end", "gene", "log2"]
        )

    def _unsorted_bins(self):
        """Descending starts on one chromosome: what both guards must catch."""
        return self._bins(
            [
                ["chr1", 300, 400, "D", 0.0],
                ["chr1", 0, 100, "A", 0.0],
                ["chr1", 100, 200, "B", 0.0],
            ]
        )

    def test_translate_region_to_bins_rejects_unsorted(self):
        bins = self._unsorted_bins()
        with self.assertRaises(ValueError) as cm:
            plots.translate_region_to_bins("chr1:0-250", bins)
        self.assertIn("sorted by start position", str(cm.exception))

    def test_update_binwise_positions_rejects_unsorted(self):
        bins = self._unsorted_bins()
        segs = self._bins([["chr1", 0, 400, "-", 0.0]])
        with self.assertRaises(ValueError):
            plots.update_binwise_positions(bins, segs)

    def test_chromosome_major_array_is_not_unsorted(self):
        """Coordinates restart at each chromosome, so a whole-genome array is
        non-monotonic by construction. The guard is per chromosome; a whole-
        table check would reject every normal input."""
        bins = self._bins(
            [
                ["chr1", 100, 200, "A", 0.0],
                ["chr1", 200, 300, "B", 0.0],
                ["chr2", 0, 100, "C", 0.0],
                ["chr2", 100, 200, "D", 0.0],
            ]
        )
        self.assertFalse(bins.data.start.is_monotonic_increasing)
        # Both chr2 bins start below 150, and are numbered from 0 on their own
        # chromosome -- chr1's two rows do not shift them.
        self.assertEqual(
            plots.translate_region_to_bins("chr2:0-150", bins),
            ("chr2", 0, 2),
        )
        cnarr, _s, _v, _x = plots.update_binwise_positions(bins)
        self.assertEqual(cnarr.start.tolist(), [0, 1, 0, 1])

    def test_absent_chromosome_yields_empty_span(self):
        bins = self._bins([["chr1", 100, 200, "A", 0.0]])
        self.assertEqual(
            plots.translate_region_to_bins("chrZZ:1-100", bins), ("chrZZ", 0, 0)
        )

    def test_sorted_ordinals_are_unchanged(self):
        """Pin the ordinals themselves, not just the absence of a raise.

        Both endpoints are a 'left' search of the *start* column. Queries that
        land exactly on a bin start cannot tell that apart from an overlap
        query -- all three agree there -- so this uses coordinates strictly
        inside a bin, where `skgenome.intersect`'s modes part company:
        'inner' would give (51, 150) and 'outer' (50, 151). Neither reproduces
        the pair, so a reroute through the overlap machinery fails here.
        """
        cnarr = cnvlib.read("formats/amplicon.cnr")
        segarr = cnvlib.read("formats/amplicon.cns")
        # Midway into chr19's 50th and 150th bins, of 273
        self.assertEqual(
            plots.translate_region_to_bins("chr19:9007547-9067688", cnarr),
            ("chr19", 51, 151),
        )
        # On a bin start, the ordinal is that bin's own index
        self.assertEqual(
            plots.translate_region_to_bins("chr19:9007428-9067555", cnarr),
            ("chr19", 50, 150),
        )
        _c, segs, _v, _x = plots.update_binwise_positions(cnarr, segarr)
        chr19 = segs.data.loc[segs.chromosome == "chr19", ["start", "end"]]
        self.assertEqual(chr19.to_numpy().tolist(), [[0, 7], [7, 25], [25, 273]])


class ByBinVariantTests(unittest.TestCase):
    """Placing variants on the bin axis.

    `scatter --by-bin -v` crashed for the whole life of the feature, so there
    is no prior rendering to preserve; these pin the intended one.
    """

    @staticmethod
    def _arrays():
        # Bins 0 and 1 abut; bin 2 is separated from them by a 300-400 gap.
        bins = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 200, "A", 0.0],
                ["chr1", 200, 300, "B", 0.0],
                ["chr1", 400, 500, "C", 0.0],
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        # Two co-binned, one in the gap, one on-target, one off each end.
        snvs = vary.VariantArray.from_rows(
            [
                ["chr1", 50, 51, "A", "G", 0.5],  # before the first bin
                ["chr1", 150, 151, "C", "T", 0.5],  # inside bin 0
                ["chr1", 160, 161, "G", "A", 0.5],  # inside bin 0, sharing it
                ["chr1", 350, 351, "T", "C", 0.5],  # in the gap after bin 1
                ["chr1", 450, 451, "A", "T", 0.5],  # inside bin 2
                ["chr1", 600, 601, "C", "G", 0.5],  # past the last bin
            ],
            columns=["chromosome", "start", "end", "ref", "alt", "alt_freq"],
        )
        return bins, snvs

    def test_variants_take_their_own_bin_and_gaps_take_the_boundary(self):
        """A position inside a bin gets that bin, not its rank among starts.

        searchsorted's default 'left' rank is one too high for anything not
        landing exactly on a bin start, which is nearly every variant. A
        variant between two bins belongs to neither, and the gap has no width
        on this axis, so it goes to the boundary with the following bin --
        index 2 below -- rather than to the preceding one, which a plain
        ``searchsorted(side="right") - 1`` would hand it.
        """
        bins, snvs = self._arrays()
        _c, _s, out, _x = plots.update_binwise_positions(bins, None, snvs)
        self.assertEqual(out.start.tolist(), [0, 0, 2, 2])

    def test_off_panel_variants_are_dropped(self):
        bins, snvs = self._arrays()
        _c, _s, out, _x = plots.update_binwise_positions(bins, None, snvs)
        self.assertEqual(len(out), 4)
        # The two with nowhere to sit, not two arbitrary rows: ref alone does
        # not distinguish them, since dropping the last two also leaves ACGT.
        self.assertEqual(
            list(zip(out.data.ref, out.data.alt, strict=True)),
            [("C", "T"), ("G", "A"), ("T", "C"), ("A", "T")],
        )

    def test_start_stays_an_integer_ordinal(self):
        """The offset is a rendering artifact and must not enter a coordinate.

        GenomicArray recasts its coordinate columns to int on every rewrap, so
        a fraction stored in `start` is destroyed by the next `by_chromosome`.
        """
        bins, snvs = self._arrays()
        _c, _s, out, _x = plots.update_binwise_positions(bins, None, snvs)
        self.assertTrue(np.issubdtype(out.data.start.dtype, np.integer))
        self.assertIn("bin_offset", out)

    def test_offsets_survive_a_rewrap_and_centre_on_the_bin(self):
        """A lone variant lands on the bin midpoint, where its log2 dot is.

        The VAF panel shares an x-axis with the log2 panel, which draws each
        bin at 0.5 * (start + end). Offsets of j/n would put every lone
        variant half a bin to the left of its own coverage point.
        """
        bins, snvs = self._arrays()
        _c, _s, out, _x = plots.update_binwise_positions(bins, None, snvs)
        per_chrom = dict(out.by_chromosome())["chr1"]
        x = plots.binwise_x(per_chrom)
        self.assertEqual(sorted(x.tolist()), [0.25, 0.75, 2.25, 2.75])

    def test_offsets_never_cross_a_bin_boundary(self):
        """The offset never carries a variant out of the bin it belongs to.

        `binwise_x` draws the dot at start plus the offset, so an offset
        reaching 1.0 would put it under the next bin's log2 point and break
        the alignment. Uses a crowded bin the other tests do not cover.
        """
        bins = cnary.CopyNumArray.from_rows(
            [["chr1", 0, 100, "A", 0.0]],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        snvs = vary.VariantArray.from_rows(
            [["chr1", 10 * i, 10 * i + 1, "A", "G", 0.5] for i in range(1, 8)],
            columns=["chromosome", "start", "end", "ref", "alt", "alt_freq"],
        )
        _c, _s, out, _x = plots.update_binwise_positions(bins, None, snvs)
        offsets = out["bin_offset"].to_numpy()
        self.assertEqual(len(offsets), 7)
        self.assertTrue(((offsets >= 0) & (offsets < 1)).all())
        self.assertAlmostEqual(offsets.mean(), 0.5)

    def test_point_cloud_does_not_depend_on_row_order(self):
        """Co-binned variants are grouped by value, not by consecutive run.

        `itertools.groupby` only sees runs, so this interleaves the two pairs
        that share a bin -- 150, 350, 160, 450. Reversing the rows would not
        do: it keeps each pair adjacent, and run-grouping would still find
        them and pass.
        """
        bins, snvs = self._arrays()
        interleaved = snvs.as_dataframe(
            snvs.data.iloc[[0, 1, 3, 2, 4, 5]].reset_index(drop=True)
        )
        _c, _s, a, _x = plots.update_binwise_positions(bins, None, snvs)
        _c, _s, b, _x = plots.update_binwise_positions(bins, None, interleaved)
        x_a, x_b = plots.binwise_x(a), plots.binwise_x(b)
        self.assertEqual(sorted(x_a.tolist()), sorted(x_b.tolist()))
        # Run-grouping would leave both pairs stacked, i.e. only 2 distinct x
        self.assertEqual(len(set(x_b.tolist())), 4)

    def test_scatter_by_bin_with_variants_renders(self):
        """Regression for the crash: int64 += float64, then two more layers."""
        # lazy: defer matplotlib import to keep headless test collection fast
        from matplotlib import pyplot  # noqa: PLC0415

        cnarr = cnvlib.read("formats/p2-20_1.cnr")
        variants = cmdutil.load_het_snps(
            "formats/na12878_na12882_mix.vcf", None, None, 0, None
        )
        fig = scatter.do_scatter(cnarr, variants=variants, by_bin=True)
        self.assertIsNotNone(fig)
        pyplot.close(fig)


class ByBinSegmentMembershipTests(unittest.TestCase):
    """Which segment a variant's B-allele frequency belongs to, under --by-bin."""

    @staticmethod
    def _arrays():
        # Bins 0 and 1 abut and carry segment A; bins 2 and 3 abut and carry
        # segment B; 300-400 is a gap in coverage that neither segment spans.
        bins = cnary.CopyNumArray.from_rows(
            [
                ["chr1", 100, 200, "A", 0.0],
                ["chr1", 200, 300, "A", 0.0],
                ["chr1", 400, 500, "B", 0.0],
                ["chr1", 500, 600, "B", 0.0],
            ],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        # Segment A overhangs the first bin, as a .cns from another tool may.
        segs = cnary.CopyNumArray.from_rows(
            [["chr1", 50, 300, "segA", 0.0], ["chr1", 400, 600, "segB", 0.0]],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        # Distinct frequencies per group, so a mis-grouping moves a median
        # rather than merely reordering equal values.
        snvs = vary.VariantArray.from_rows(
            [
                ["chr1", 60, 61, "A", "G", 0.45],  # before bin 0, inside segment A
                ["chr1", 150, 151, "A", "G", 0.20],  # bin 0, segment A
                ["chr1", 250, 251, "A", "G", 0.40],  # bin 1, segment A
                ["chr1", 320, 321, "A", "G", 0.40],  # in the gap, no segment
                ["chr1", 340, 341, "A", "G", 0.40],  # in the gap, no segment
                ["chr1", 450, 451, "A", "G", 0.10],  # bin 2, segment B
                ["chr1", 550, 551, "A", "G", 0.10],  # bin 3, segment B
            ],
            columns=["chromosome", "start", "end", "ref", "alt", "alt_freq"],
        )
        return bins, segs, snvs

    @staticmethod
    def _levels(variants, segments):
        return sorted(
            (seg.gene, round(float(value), 12))
            for seg, value in scatter.get_segment_vafs(variants, segments)
        )

    def test_gap_variants_join_no_segment(self):
        """The gap pair sits on bin 2, inside segment B's span on that axis.

        Genomically it is in neither segment, so segment B's level must stay
        the median of its own two variants. Grouping on the ordinal would
        pull the 0.40 pair in and pull the level to 0.25.
        """
        bins, segs, snvs = self._arrays()
        _c, out_segs, out_snvs, _x = plots.update_binwise_positions(bins, segs, snvs)
        self.assertEqual(out_snvs.start.tolist(), [0, 1, 2, 2, 2, 3])
        self.assertEqual(
            self._levels(out_snvs, out_segs), [("segA", 0.30), ("segB", 0.10)]
        )

    def test_by_bin_levels_match_the_genomic_view(self):
        """Same drawn levels on either axis, over the variants both can show."""
        bins, segs, snvs = self._arrays()
        _c, out_segs, out_snvs, _x = plots.update_binwise_positions(bins, segs, snvs)
        # Semantics: the retained variants only -- a level must never
        # summarize evidence the bin axis dropped and the reader cannot see.
        retained = snvs.as_dataframe(snvs.data.loc[out_snvs.data.index])
        self.assertEqual(self._levels(out_snvs, out_segs), self._levels(retained, segs))

    def test_membership_survives_subsetting_the_segments(self):
        """Membership pairs by index label, not by row position.

        The chromosome view hands `get_segment_vafs` a windowed subset whose
        rows keep their labels but no longer sit at their original positions.
        """
        bins, segs, snvs = self._arrays()
        _c, out_segs, out_snvs, _x = plots.update_binwise_positions(bins, segs, snvs)
        later = out_segs[out_segs.start >= 2]
        self.assertEqual(later.data.index.tolist(), [1])
        self.assertEqual(self._levels(out_snvs, later), [("segB", 0.10)])

    def test_unusable_segment_labels_are_refused(self):
        """The join needs unique labels, and -1 is the no-segment fill.

        Neither arises from a file, but an API caller can build one, and the
        failure would otherwise be a silently merged or poisoned level.
        """
        bins, segs, snvs = self._arrays()
        for labels in ([0, 0], [-1, 1]):
            bad = segs.copy()
            bad.data.index = pd.Index(labels)
            with self.assertRaises(ValueError) as cm:
                plots.update_binwise_positions(bins, bad, snvs)
            self.assertIn("unique index labels", str(cm.exception))


class SegmentVafNaNTests(unittest.TestCase):
    """A NaN allele frequency belongs to neither side of the 0.5 split.

    skgenome/tabio/vcfio.py leaves alt_freq as NaN when a heterozygous call
    has no usable depth or allele counts: the fraction is unknown, not zero
    (#407).
    """

    @staticmethod
    def _rows(freqs):
        return vary.VariantArray.from_rows(
            [
                ["chr1", 100 + 10 * i, 101 + 10 * i, "A", "G", freq]
                for i, freq in enumerate(freqs)
            ],
            columns=["chromosome", "start", "end", "ref", "alt", "alt_freq"],
        )

    @staticmethod
    def _segment():
        return cnary.CopyNumArray.from_rows(
            [["chr1", 0, 1000, "segA", 0.0]],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )

    def _levels(self, freqs):
        return sorted(
            round(float(value), 12)
            for _seg, value in scatter.get_segment_vafs(
                self._rows(freqs), self._segment()
            )
        )

    def test_nan_does_not_poison_the_lower_level(self):
        """One unknown frequency erased the whole below-0.5 trend line.

        np.median of a group holding any NaN is NaN, and matplotlib draws
        nothing for it, so the segment silently lost a line.
        """
        without_nan = self._levels([0.2, 0.3, 0.4, 0.7, 0.8])
        self.assertEqual(without_nan, [0.3, 0.75])
        # The NaN row must change neither level, on either side.
        self.assertEqual(
            self._levels([0.2, 0.3, 0.4, float("nan"), 0.7, 0.8]), without_nan
        )

    def test_nan_does_not_count_toward_the_two_member_minimum(self):
        """Discriminates the other plausible fix, np.nanmedian.

        Taking the median with NaN skipped but the group unchanged would
        draw a level here from a single variant, which the two-member
        minimum exists to refuse.
        """
        self.assertEqual(self._levels([0.2, float("nan"), 0.7, 0.8]), [0.75])

    def test_an_all_nan_group_draws_nothing(self):
        self.assertEqual(self._levels([float("nan")] * 4), [])
