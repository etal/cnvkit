#!/usr/bin/env python
"""Unit tests for the CNVkit library, cnvlib."""
from __future__ import absolute_import, division, print_function

import unittest

from skgenome import tabio


class IOTests(unittest.TestCase):
    """Tests for I/O modules."""

    def test_empty(self):
        """Instantiate from an empty file."""
        for fmt in ("auto", "bed", #"bed3", "bed4",
                    "interval", "tab", "text"):
            regions = tabio.read("formats/empty", fmt=fmt)
            self.assertEqual(len(regions), 0)

    def test_read_auto(self):
        for fname, nrows in (("formats/empty", 0),
                             ("formats/amplicon.bed", 1433),
                             ("formats/amplicon.text", 1433),
                             ("formats/nv2_baits.interval_list", 6809),
                             ("formats/refflat-mini.txt", 100),
                             ("formats/example.gff", 6),
                            ):
            self.assertEqual(len(tabio.read_auto(fname)), nrows)
            with open(fname) as handle:
                self.assertEqual(len(tabio.read_auto(handle)), nrows)

    def test_read_bed(self):
        """Read the BED format."""
        fname = "formats/amplicon.bed"
        regions = tabio.read(fname, "bed")
        self.assertEqual(len(regions), linecount(fname))
        self.assertEqual(regions.sample_id, "amplicon")

    def test_read_gff(self):
        """Read the GFF format."""
        for fname, nrows, sample_id in (
            ("formats/example.gff", 6, "example"),
            ("formats/GRCh37_BRAF.gff.gz", 49, "GRCh37_BRAF"),
        ):
            regions = tabio.read(fname, 'gff')
            self.assertEqual(len(regions), nrows)
            self.assertEqual(regions.sample_id, sample_id)

    def test_read_ilist(self):
        """Read the interval list format."""
        regions = tabio.read("formats/nv2_baits.interval_list", "interval")
        self.assertEqual(len(regions), 6809)
        self.assertEqual(regions.sample_id, "nv2_baits")

    def test_read_picardhs(self):
        """Read Picard CalculateHsMetrics PER_TARGET_COVERAGE format."""
        fname = "picard/p2-5_5.antitargetcoverage.csv"
        cna = tabio.read(fname, "picardhs")
        self.assertEqual(len(cna), linecount(fname) - 1)
        self.assertEqual(cna.sample_id, "p2-5_5")

    def test_read_refflat(self):
        """Read the UCSC 'refFlat' format."""
        fname = "formats/refflat-mini.txt"
        regions = tabio.read(fname, 'refflat')
        self.assertEqual(len(regions), linecount(fname))
        self.assertEqual(13, regions.chromosome.nunique())

    def test_read_seg(self):
        """Read the SEG format."""
        for fname, header_len, args in (
            # Convert integer chrom. IDs to hg19 names
            ('formats/cw-tr-log2.seg', 1,
             ({'23': 'X', '24': 'Y', '25': 'M'}, "chr", False)),
            # Convert segmented array CGH data in log10 scale to log2
            ('formats/acgh-log10.seg', 1, (None, None, True)),
            # From DNAcopy, with a stray warning message from R
            ('formats/warning.seg', 2, (None, None, False)),
        ):
            expect_lines = linecount(fname) - header_len
            seen_lines = 0
            for _sample_id, dframe in tabio.seg.parse_seg(fname, *args):
                seen_lines += len(dframe)
            self.assertEqual(seen_lines, expect_lines)

    def test_read_text(self):
        """Read the text region format."""
        fname = "formats/amplicon.text"
        regions = tabio.read(fname, "text")
        self.assertEqual(len(regions), linecount(fname))
        self.assertEqual(regions.sample_id, "amplicon")

    def test_read_vcf(self):
        """Read the VCF format."""
        # Paired VCF with full info
        fname = "formats/na12878_na12882_mix.vcf"
        v1 = tabio.read(fname, "vcf")
        self.assertLess(len(v1), linecount(fname))
        self.assertLess(0, len(v1))
        for sid in ("NA12882", "NA12878"):
            v2 = tabio.read(fname, "vcf", sample_id=sid)
            self.assertEqual(v2.sample_id, sid)
            self.assertEqual(len(v1), len(v2))
        for kwarg in ({'min_depth': 100},
                      {'skip_somatic': True},
                      {'skip_reject': True}):
            v3 = tabio.read(fname, "vcf", **kwarg)
            self.assertLess(len(v3), len(v1))
            self.assertLess(0, len(v3),
                            "%d variants left after filter %r"
                            % (len(v3), list(kwarg)[0]))
        # VCF header, no samples, no records
        v4 = tabio.read('formats/nosample.vcf', 'vcf')
        self.assertEqual(len(v4), 0)
        self.assertEqual(v4.sample_id, 'nosample')
        # VCF with 1 sample, no records
        v5 = tabio.read('formats/blank.vcf', 'vcf', sample_id='Blank')
        self.assertEqual(len(v5), 0)
        self.assertEqual(v5.sample_id, 'Blank')


# == helpers ==

def linecount(filename):
    i = -1
    with open(filename) as handle:
        for i, _line in enumerate(handle):
            pass
        return i + 1


if __name__ == '__main__':
    unittest.main(verbosity=2)
