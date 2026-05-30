#!/usr/bin/env python
"""Tests for cnvlib.samutil, focused on BAM index (BAI vs. CSI) selection."""

import logging
import os
import tempfile
import unittest

import pysam

from cnvlib import samutil
from cnvlib.samutil import _BAI_MAX_CONTIG_LENGTH as BAI_LIMIT

logging.basicConfig(level=logging.ERROR, format="%(message)s")


def write_sorted_bam(path, contig_length, read_start):
    """Write a coordinate-sorted single-read BAM with one contig.

    The contig is declared `contig_length` bp long (only the header carries the
    length; no reference sequence is stored), and one mapped read is placed at
    `read_start`. Placing a read near the end of a >2**29 contig is what makes a
    BAI index actually impossible, so this exercises the real failure mode.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrBig", "LN": contig_length}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        read = pysam.AlignedSegment(out.header)
        read.query_name = "read1"
        read.query_sequence = "ACGT" * 10
        read.flag = 0
        read.reference_id = 0
        read.reference_start = read_start
        read.mapping_quality = 60
        read.cigarstring = "40M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 40)
        out.write(read)


class IndexSelectionTests(unittest.TestCase):
    """ensure_bam_index picks BAI for short contigs and CSI for long ones."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.bam = os.path.join(self.tmpdir, "sample.bam")

    def tearDown(self):
        for name in os.listdir(self.tmpdir):
            os.remove(os.path.join(self.tmpdir, name))
        os.rmdir(self.tmpdir)

    def test_short_contig_builds_bai(self):
        """A normal (sub-2**29) contig is indexed as BAI, unchanged behavior."""
        write_sorted_bam(self.bam, contig_length=100_000, read_start=1000)
        self.assertFalse(samutil._bam_needs_csi(self.bam, pysam))
        index_fname = samutil.ensure_bam_index(self.bam)
        self.assertEqual(index_fname, self.bam + ".bai")
        self.assertTrue(os.path.isfile(index_fname))
        self.assertFalse(os.path.isfile(self.bam + ".csi"))

    def test_long_contig_builds_usable_csi(self):
        """A >2**29 contig is auto-indexed as CSI and remains queryable."""
        write_sorted_bam(self.bam, contig_length=600_000_000, read_start=599_000_000)
        self.assertTrue(samutil._bam_needs_csi(self.bam, pysam))
        index_fname = samutil.ensure_bam_index(self.bam)
        self.assertEqual(index_fname, self.bam + ".csi")
        self.assertTrue(os.path.isfile(index_fname))
        self.assertFalse(os.path.isfile(self.bam + ".bai"))
        # The CSI index must actually work for downstream pysam operations.
        with pysam.AlignmentFile(self.bam, "rb") as af:
            self.assertEqual(sum(1 for _ in af.fetch("chrBig")), 1)
        idxstats = pysam.idxstats(self.bam)
        self.assertIn("chrBig\t600000000\t1\t0", idxstats)

    def test_contig_exactly_at_limit_uses_bai(self):
        """A contig of length exactly 2**29 still fits BAI (strict-> threshold)."""
        write_sorted_bam(self.bam, contig_length=BAI_LIMIT, read_start=BAI_LIMIT - 1000)
        self.assertFalse(samutil._bam_needs_csi(self.bam, pysam))
        index_fname = samutil.ensure_bam_index(self.bam)
        self.assertEqual(index_fname, self.bam + ".bai")
        self.assertTrue(os.path.isfile(index_fname))

    def test_existing_fresh_csi_is_reused(self):
        """A CSI index newer than the BAM is reused without rebuilding."""
        write_sorted_bam(self.bam, contig_length=600_000_000, read_start=599_000_000)
        pysam.index("-c", self.bam)
        csi = self.bam + ".csi"
        # The freshly built CSI is newer than the BAM, so it must be reused as-is;
        # an unwanted rebuild would rewrite the file and change its mtime.
        before_ns = os.stat(csi).st_mtime_ns
        index_fname = samutil.ensure_bam_index(self.bam)
        self.assertEqual(index_fname, csi)
        self.assertEqual(os.stat(csi).st_mtime_ns, before_ns)

    def test_stale_bai_removed_when_csi_built(self):
        """A leftover .bai is deleted after a CSI index is built for the BAM."""
        write_sorted_bam(self.bam, contig_length=600_000_000, read_start=599_000_000)
        # A leftover .bai from a prior alignment; the CSI rebuild must remove it.
        stale_bai = self.bam + ".bai"
        with open(stale_bai, "wb") as handle:
            handle.write(b"BAI\x01")
        index_fname = samutil.ensure_bam_index(self.bam)
        self.assertEqual(index_fname, self.bam + ".csi")
        self.assertFalse(os.path.isfile(stale_bai))


if __name__ == "__main__":
    unittest.main(verbosity=2)
