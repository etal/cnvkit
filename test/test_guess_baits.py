#!/usr/bin/env python
"""Tests for the guess_baits.py depth-scan parser core.

These tests are samtools-free: ``subprocess.Popen`` is mocked so the parser
iterates synthetic ``samtools depth`` text lines and never touches a BAM.
``scan_depth`` now drives ``Popen`` as a context manager and checks its
``returncode``, so the mock wires stdout onto the ``__enter__`` target and
supplies a ``returncode``.

Includes a regression guard for #542, where a bytes-mode split delimiter
(``line.split(b"\\t")`` over a bytes ``proc.stdout``) silently produced bytes
chromosome names (``b'chr1'``) in the emitted regions, which then leaked into
the written BED/tab output.
"""

import unittest
from unittest import mock

import numpy as np
import pandas as pd
from pandas.api.types import is_float_dtype, is_integer_dtype, is_object_dtype

from cnvlib.cli import guess_baits
from cnvlib.cli.guess_baits import Region


def scan_lines(lines, bam_fnames, min_depth, returncode=0):
    """Run scan_depth over synthetic ``samtools depth`` text lines.

    ``cnvlib.cli.guess_baits.subprocess.Popen`` is replaced and used as a
    context manager: ``proc = popen.return_value.__enter__.return_value`` gets
    a str-line ``stdout`` iterator (as the parser expects after the
    ``text=True`` fix) and a ``returncode`` (0 = success). The parser iterates
    whatever the mocked stdout yields, so these synthetic str lines are
    precisely what the code under test parses.
    """
    with mock.patch("cnvlib.cli.guess_baits.subprocess.Popen") as popen:
        proc = popen.return_value.__enter__.return_value
        proc.stdout = iter(lines)
        proc.returncode = returncode
        return list(guess_baits.scan_depth("access.bed", bam_fnames, min_depth))


class ScanDepthTests(unittest.TestCase):
    """scan_depth: parsing, coordinates, edge-trimming, multi-sample scaling."""

    def test_chromosome_is_str_regression_542(self):
        """#542: emitted chromosome must be str, not bytes, and Popen text=True.

        On the buggy ``line.split(b"\\t")`` over str stdout this raises
        TypeError; on the fixed ``line.split("\\t")`` it yields "chr1". The
        parser only ever sees str lines because ``Popen(..., text=True)``
        decodes samtools stdout; dropping ``text=True`` reverts real output to
        bytes, so this test also pins that constructor argument. The mock is
        inlined here (rather than via scan_lines) to inspect popen.call_args.
        """
        # Trailing low-depth line flushes the final enriched run.
        lines = ["chr1\t1000\t7\n", "chr1\t1001\t0\n"]
        with mock.patch("cnvlib.cli.guess_baits.subprocess.Popen") as popen:
            proc = popen.return_value.__enter__.return_value
            proc.stdout = iter(lines)
            proc.returncode = 0
            regions = list(guess_baits.scan_depth("access.bed", ["a.bam"], 1))
        self.assertEqual(len(regions), 1)
        chrom = regions[0].chromosome
        self.assertIsInstance(chrom, str)
        self.assertEqual(chrom, "chr1")
        self.assertNotEqual(chrom, b"chr1")
        # text=True makes samtools stdout str lines, not bytes.
        self.assertTrue(popen.call_args.kwargs.get("text"))

    def test_coordinates_and_edge_trimming(self):
        """1-indexed samtools pos -> 0-indexed half-open; trim to >= 0.5*peak.

        Depth profile at samtools pos 1000..1003 is [4, 8, 12, 4]. All four
        bases clear min_depth=3 and form one run, but peak=12 so the half-peak
        cutoff is 6: only the [8, 12] bases survive trimming. The reported
        depth is the mean over the *trimmed* window (10.0), not the peak (12),
        which distinguishes mean-over-window from peak reporting.
        """
        lines = [
            "chr1\t1000\t4\n",
            "chr1\t1001\t8\n",
            "chr1\t1002\t12\n",
            "chr1\t1003\t4\n",
            "chr1\t1004\t0\n",  # flush
        ]
        regions = scan_lines(lines, ["a.bam"], 3)
        self.assertEqual(len(regions), 1)
        reg = regions[0]
        # samtools pos 1001 (depth 8) -> 0-indexed start 1000;
        # last kept base is pos 1002 (depth 12) -> half-open end 1002.
        self.assertEqual(reg.start, 1000)
        self.assertEqual(reg.end, 1002)
        self.assertEqual(reg.depth, 10.0)

    def test_multisample_sum_and_threshold_scaling(self):
        """Two BAMs: depths summed across trailing columns; min_depth *= nsamples.

        Two samples scale min_depth=3 to an effective threshold of 6. Position
        1000 sums to 5 -- it clears the unscaled threshold 3 but fails the
        scaled 6, so it never joins the run; position 1001 sums to 8 and is the
        only enriched base. The result is exactly one
        Region("chr5", 1000, 1001, 8.0).

        This bites the scaling contract: dropping ``min_depth *= len(bam_fnames)``
        in guess_baits.scan_depth would admit position 1000 as well, flipping
        the region to Region("chr5", 999, 1001, 6.5) -- both the start
        (999 vs 1000) and the depth (6.5 vs 8.0) change. Depth 8.0 (= 4 + 4)
        also proves the per-base sum across the two trailing samtools columns.
        """
        lines = [
            "chr5\t1000\t2\t3\n",  # sum 5 < scaled 6 (but > unscaled 3)
            "chr5\t1001\t4\t4\n",  # sum 8 >= scaled 6
            "chr5\t1002\t0\t0\n",  # flush
        ]
        regions = scan_lines(lines, ["a.bam", "b.bam"], 3)
        self.assertEqual(len(regions), 1)
        reg = regions[0]
        # Only pos 1001 (sum 8) enriched: start 1000, half-open end 1001.
        self.assertEqual(reg.start, 1000)
        self.assertEqual(reg.end, 1001)
        self.assertEqual(reg.depth, 8.0)

    def test_chromosome_change_flushes_region(self):
        """A chromosome change flushes the current enriched run."""
        lines = [
            "chr1\t100\t9\n",
            "chr2\t100\t0\n",  # new chrom flushes the chr1 run
        ]
        regions = scan_lines(lines, ["a.bam"], 5)
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].chromosome, "chr1")
        self.assertEqual(regions[0].start, 99)
        self.assertEqual(regions[0].end, 100)

    def test_final_run_flushed_at_eof(self):
        """A run enriched through the LAST line is flushed at end-of-stream.

        There is no trailing sub-threshold line to close the run, so the region
        exists only if end-of-stream flushes it. The pre-fix parser dropped it,
        yielding 0 regions.
        """
        lines = [
            "chr1\t500\t10\n",
            "chr1\t501\t10\n",  # enriched through EOF, no closing line
        ]
        regions = scan_lines(lines, ["a.bam"], 5)
        self.assertEqual(len(regions), 1)
        reg = regions[0]
        self.assertEqual(reg.chromosome, "chr1")
        self.assertEqual(reg.start, 499)
        self.assertEqual(reg.end, 501)
        self.assertEqual(reg.depth, 10.0)

    def test_chromosome_change_keeps_first_base_of_new_run(self):
        """A back-to-back enriched chrom change keeps the new run's first base.

        chr1 and chr2 are both enriched with no intervening sub-threshold line;
        a third chrom flushes chr2. The chr2 run must start at its own first
        position (pos-1), not one base in -- the pre-fix parser dropped the
        first base after a chromosome change.
        """
        lines = [
            "chr1\t100\t10\n",
            "chr2\t200\t10\n",  # enriched, immediately after chr1
            "chr3\t300\t0\n",  # flush chr2
        ]
        regions = scan_lines(lines, ["a.bam"], 5)
        self.assertEqual(len(regions), 2)
        self.assertEqual(regions[0].chromosome, "chr1")
        self.assertEqual((regions[0].start, regions[0].end), (99, 100))
        chr2 = regions[1]
        self.assertEqual(chr2.chromosome, "chr2")
        # First base of the chr2 run is kept: samtools pos 200 -> start 199.
        self.assertEqual(chr2.start, 199)
        self.assertEqual(chr2.end, 200)

    def test_depth_zero_hole_breaks_run(self):
        """A gap in reported positions (an omitted depth-0 base) breaks the run.

        samtools omits zero-depth bases, so pos 101 is absent between the
        enriched pos 100 and pos 102. Positions must be consecutive to extend a
        run, so this must produce two separate 1-base regions [99,100) and
        [101,102) -- not one drifted [99,101) that silently spans the hole.
        """
        lines = [
            "chr1\t100\t10\n",
            "chr1\t102\t10\n",  # pos 101 omitted (depth 0): gap breaks the run
            "chr1\t200\t0\n",  # flush
        ]
        regions = scan_lines(lines, ["a.bam"], 5)
        self.assertEqual(len(regions), 2)
        self.assertEqual((regions[0].start, regions[0].end), (99, 100))
        self.assertEqual((regions[1].start, regions[1].end), (101, 102))

    def test_nonzero_returncode_raises(self):
        """A nonzero samtools exit code after the stream raises RuntimeError."""
        lines = ["chr1\t100\t0\n"]  # nothing enriched; failure is the contract
        with self.assertRaisesRegex(RuntimeError, r"exit 1.*access\.bed"):
            scan_lines(lines, ["a.bam"], 5, returncode=1)

    def test_scan_depth_empty_input_returns_empty_frame(self):
        """_scan_depth over an all-sub-threshold stream yields an empty frame.

        The full drop_small(merge_gaps(scan_depth(...))) pipeline must survive
        an empty region list and return a DataFrame with 0 rows and the Region
        field columns -- no StopIteration/RuntimeError from empty generators.
        """
        lines = ["chr1\t100\t1\n", "chr1\t101\t2\n"]  # all below threshold
        with mock.patch("cnvlib.cli.guess_baits.subprocess.Popen") as popen:
            proc = popen.return_value.__enter__.return_value
            proc.stdout = iter(lines)
            proc.returncode = 0
            bed_fname, frame = guess_baits._scan_depth(("bed", ["a.bam"], 5, 25, 10))
        self.assertEqual(bed_fname, "bed")
        self.assertEqual(len(frame), 0)
        self.assertEqual(list(frame.columns), list(Region._fields))
        # merge_gaps is empty-safe on its own, too.
        self.assertEqual(list(guess_baits.merge_gaps(iter([]), 25)), [])

    def test_empty_and_nonempty_concat_keeps_numeric_dtypes(self):
        """Mixed empty+non-empty _scan_depth frames concat without object dtype.

        ``pd.DataFrame.from_records([], columns=...)`` infers all-object dtype
        for an empty chunk, and scan_targets concats every chunk's frame. If an
        empty frame's columns stayed object, pd.concat would down-cast the
        numeric columns to object across the whole result, and the
        ``guess-baits --coverage`` path's ``np.log2(depth / depth.median())``
        would raise TypeError. Genome-wide ``--access`` yields mostly-empty
        chunks, so this mixed concat is the common case, not a corner case.

        _scan_depth now coerces the empty frame to Region's dtypes, so the
        concatenated start/end stay integer, depth stays float (not object),
        and the log2 normalization runs. This fails on the pre-fix code.
        """

        def scan_depth_frame(lines, min_depth, min_gap, min_length):
            with mock.patch("cnvlib.cli.guess_baits.subprocess.Popen") as popen:
                proc = popen.return_value.__enter__.return_value
                proc.stdout = iter(lines)
                proc.returncode = 0
                _, frame = guess_baits._scan_depth(
                    ("bed", ["a.bam"], min_depth, min_gap, min_length)
                )
            return frame

        # All sub-threshold -> empty chunk; one enriched run -> non-empty chunk.
        empty = scan_depth_frame(["chr1\t100\t1\n", "chr1\t101\t2\n"], 5, 25, 10)
        nonempty = scan_depth_frame(
            ["chr1\t1000\t10\n", "chr1\t1001\t10\n", "chr1\t1002\t0\n"], 1, 25, 1
        )
        self.assertEqual(len(empty), 0)
        self.assertEqual(len(nonempty), 1)

        cat = pd.concat([empty, nonempty])
        self.assertTrue(is_integer_dtype(cat["start"]))
        self.assertTrue(is_integer_dtype(cat["end"]))
        self.assertTrue(is_float_dtype(cat["depth"]))
        self.assertFalse(is_object_dtype(cat["depth"]))
        # The actual --coverage failure mode: object dtype breaks this ufunc.
        ratio = np.log2(cat["depth"] / cat["depth"].median())
        self.assertTrue(np.isfinite(ratio).all())


class MergeGapsTests(unittest.TestCase):
    """merge_gaps: strict ``<`` gap merge, chromosome-aware, residual flush."""

    def test_small_gap_merged_large_gap_kept(self):
        regions = [
            Region("c", 100, 200, 1.0),
            Region("c", 250, 300, 1.0),  # gap 50 from prev
            Region("c", 500, 600, 1.0),  # gap 200 from prev
        ]
        merged = list(guess_baits.merge_gaps(iter(regions), 100))
        self.assertEqual(len(merged), 2)
        # First two merged: end extended to 300; last kept separate (residual).
        self.assertEqual((merged[0].start, merged[0].end), (100, 300))
        self.assertEqual((merged[1].start, merged[1].end), (500, 600))

    def test_gap_equal_to_min_gap_not_merged(self):
        """Boundary: gap == min_gap is NOT merged (strict ``<``)."""
        regions = [
            Region("c", 100, 200, 1.0),
            Region("c", 250, 300, 1.0),  # gap exactly 50
        ]
        merged = list(guess_baits.merge_gaps(iter(regions), 50))
        self.assertEqual(len(merged), 2)
        self.assertEqual((merged[0].start, merged[0].end), (100, 200))
        self.assertEqual((merged[1].start, merged[1].end), (250, 300))

    def test_different_chromosomes_never_merged(self):
        """Regions on different chromosomes are kept separate.

        The coordinate difference here (5 - 1000 = -995) would satisfy a naive
        ``diff < min_gap`` test; without the chromosome guard the pre-fix code
        merged them into one corrupt region with end < start.
        """
        regions = [
            Region("chr1", 900, 1000, 1.0),
            Region("chr2", 5, 200, 1.0),
        ]
        merged = list(guess_baits.merge_gaps(iter(regions), 25))
        self.assertEqual(len(merged), 2)
        self.assertEqual(
            (merged[0].chromosome, merged[0].start, merged[0].end),
            ("chr1", 900, 1000),
        )
        self.assertEqual(
            (merged[1].chromosome, merged[1].start, merged[1].end),
            ("chr2", 5, 200),
        )
        # Neither region is corrupt (end >= start on both).
        self.assertGreaterEqual(merged[1].end, merged[1].start)

    def test_merged_depth_is_length_weighted(self):
        """Merged depth is the length-weighted mean of the sub-regions.

        (10.0 * 100 + 4.0 * 50) / (100 + 50) == 8.0, distinct from both the
        arithmetic mean (7.0) and either input depth.
        """
        regions = [
            Region("c", 100, 200, 10.0),  # length 100
            Region("c", 210, 260, 4.0),  # length 50, gap 10 < 25
        ]
        merged = list(guess_baits.merge_gaps(iter(regions), 25))
        self.assertEqual(len(merged), 1)
        self.assertEqual((merged[0].start, merged[0].end), (100, 260))
        self.assertEqual(merged[0].depth, 8.0)

    def test_single_region_passes_through_unchanged(self):
        """A lone region is emitted with its depth bit-exact (no re-averaging)."""
        regions = [Region("c", 100, 200, 7.3)]
        merged = list(guess_baits.merge_gaps(iter(regions), 25))
        self.assertEqual(len(merged), 1)
        self.assertEqual((merged[0].start, merged[0].end), (100, 200))
        self.assertEqual(merged[0].depth, 7.3)


class DropSmallTests(unittest.TestCase):
    """drop_small: keep regions with end - start >= min_length."""

    def test_filters_by_length_inclusive_boundary(self):
        regions = [
            Region("c", 100, 200, 1.0),  # length 100 -- kept (boundary)
            Region("c", 250, 300, 1.0),  # length 50  -- dropped
            Region("c", 500, 600, 1.0),  # length 100 -- kept
        ]
        kept = list(guess_baits.drop_small(iter(regions), 100))
        self.assertEqual(len(kept), 2)
        self.assertEqual((kept[0].start, kept[0].end), (100, 200))
        self.assertEqual((kept[1].start, kept[1].end), (500, 600))


if __name__ == "__main__":
    unittest.main(verbosity=2)
