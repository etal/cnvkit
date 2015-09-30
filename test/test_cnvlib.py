#!/usr/bin/env python

"""Unit tests for the CNVkit library, cnvlib."""
from __future__ import absolute_import, division, print_function

import unittest

import numpy as np
# from Bio._py3k import StringIO

import cnvlib
# Import all modules as a smoke test
from cnvlib import (antitarget, commands, core, coverage, diagram, export, fix,
                    importers, metrics, ngfrills, params, plots, reference,
                    reports, segmentation, smoothing,
                    gary, cnary, vary, rary)


class GaryTests(unittest.TestCase):

    def setUp(self):
        self.ex_cnr = cnvlib.read('formats/reference-tr.cnn')

    def test_empty(self):
        """Instantiate from an empty file."""
        garr = gary.GenomicArray.read("formats/empty")
        self.assertEqual(len(garr), 0)

    def test_iter(self):
        """Test iteration."""
        rows = iter(self.ex_cnr)
        firstrow = next(rows)
        self.assertEqual(tuple(firstrow), tuple(self.ex_cnr[0]))
        i = 0
        for i, _row in enumerate(rows):
            pass
        self.assertEqual(i + 2, len(self.ex_cnr))

    def test_copy(self):
        """Test creation of an independent copy of the object."""
        dupe = self.ex_cnr.copy()
        self.assertEqual(tuple(self.ex_cnr[3]), tuple(dupe[3]))
        self.ex_cnr[3, 'log2'] = -10.0
        self.assertNotEqual(tuple(self.ex_cnr[3]), tuple(dupe[3]))

    # def test_by_bin(self):
    def test_by_chromosome(self):
        for fname in ("formats/amplicon.cnr", "formats/cl_seq.cns"):
            cnarr = cnvlib.read(fname)
            row_count = 0
            for _chrom, rows in cnarr.by_chromosome():
                row_count += len(rows)
            self.assertEqual(row_count, len(cnarr))

    def test_select(self):
        """Test sugary selection of a subset of the data array."""
        num_bg_rows = len(self.ex_cnr[self.ex_cnr['gene'] == 'Background'])
        self.assertEqual(len(self.ex_cnr.select(gene='Background')),
                         num_bg_rows)
        selector = lambda row: row['gene'] == 'Background'
        self.assertEqual(len(self.ex_cnr.select(selector)), num_bg_rows)

    def test_shuffle_sort(self):
        """Test shuffling and re-sorting the data array."""
        orig_cvg = tuple(self.ex_cnr['log2'][:10])
        self.assertEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)
        self.ex_cnr.shuffle()
        self.assertNotEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)
        self.ex_cnr.sort()
        self.assertEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)



class CNATests(unittest.TestCase):
    """Tests for the CopyNumArray class."""

    def setUp(self):
        self.ex_cnr = cnvlib.read('formats/reference-tr.cnn')

    def test_empty(self):
        """Instantiate from an empty file."""
        cnarr = cnvlib.read("formats/empty")
        self.assertEqual(len(cnarr), 0)

    def test_basic(self):
        """Test basic container functionality and magic methods."""
        # Length
        self.assertEqual(len(self.ex_cnr), 27526)
        # Equality
        same = cnvlib.read('formats/reference-tr.cnn')
        self.assertEqual(self.ex_cnr, same)
        # Item access
        orig = self.ex_cnr[0]
        self.ex_cnr[0] = orig
        self.ex_cnr[3:4] = self.ex_cnr[3:4]
        self.ex_cnr[6:10] = self.ex_cnr[6:10]
        self.assertEqual(tuple(self.ex_cnr[0]), tuple(same[0]))
        self.assertEqual(self.ex_cnr[3:6], same[3:6])

    # def test_by_gene(self):
    # def test_by_segment(self):

    def test_center_all(self):
        """Test median-recentering."""
        chr1 = self.ex_cnr.in_range('chr1')
        self.assertAlmostEqual(0, np.median(chr1['log2']), places=1)
        chr1.center_all()
        orig_chr1_cvg = np.median(chr1['log2'])
        self.assertAlmostEqual(0, orig_chr1_cvg)
        chr1plus2 = chr1.copy()
        chr1plus2['log2'] += 2.0
        chr1plus2.center_all()
        self.assertAlmostEqual(np.median(chr1plus2['log2']), orig_chr1_cvg)

    def test_drop_extra_columns(self):
        """Test removal of optional 'gc' column."""
        self.assertTrue('gc' in self.ex_cnr)
        cleaned = self.ex_cnr.drop_extra_columns()
        self.assertTrue('gc' not in cleaned)
        self.assertTrue((cleaned['log2'] == self.ex_cnr['log2']).all())

    # def test_extend(self):
    # def test_in_range(self):

    # def test_squash_genes(self):



class RATests(unittest.TestCase):
    """Tests for RegionArray class."""

    def test_empty(self):
        """Instantiate from an empty file."""
        regions = rary.RegionArray.read("formats/empty")
        self.assertEqual(len(regions), 0)

    def test_read_bed(self):
        """Read the BED format."""
        regions = rary.RegionArray.read("formats/amplicon.bed")
        self.assertEqual(len(regions), 1433)

    def test_read_text(self):
        """Read the text region format."""
        regions = rary.RegionArray.read("formats/amplicon.text")
        self.assertEqual(len(regions), 1433)

    def test_read_ilist(self):
        """Read the interval list format."""
        regions = rary.RegionArray.read("formats/nv2_baits.interval_list")
        self.assertEqual(len(regions), 6809)



class ImporterTests(unittest.TestCase):
    """Tests for importers functionality."""

    def test_import_picard(self):
        """Test loading a Picard targetcoverage file."""
        fname = 'picard/p2-5_5.antitargetcoverage.csv'
        cna = importers.import_picard_pertargetcoverage(fname)
        self.assertTrue(len(cna) > 1)

    def test_import_seg(self):
        """Test loading SEG format."""
        for fname, args in (
            # cnvkit.py import-seg cw-tr-log2.seg -p chr -c human -d tmp
            ('formats/cw-tr-log2.seg', ({'23': 'X', '24': 'Y', '25': 'M'}, "chr", False)),
            # cnvkit.py import-seg --from-log10 acgh-log10.seg -d tmp/
            ('formats/acgh-log10.seg', (None, None, True))):
            expect_lines = linecount(fname) - 1
            seen_lines = 0
            for cns in importers.import_seg(fname, *args):
                seen_lines += len(cns)
            self.assertEqual(seen_lines, expect_lines)



class CommandTests(unittest.TestCase):
    """Tests for top-level commands."""

    def test_antitarget(self):
        """The 'antitarget' command."""
        baits_fname = "formats/nv2_baits.interval_list"
        access_fname = "../data/access-5k-mappable.hg19.bed"
        self.assertTrue(0 < len(list(
            commands.do_antitarget(baits_fname))))
        self.assertTrue(0 < len(list(
            commands.do_antitarget(baits_fname, access_fname))))
        self.assertTrue(0 < len(list(
            commands.do_antitarget(baits_fname, access_fname, 200000))))
        self.assertTrue(0 < len(list(
            commands.do_antitarget(baits_fname, access_fname, 10000, 5000))))

    def test_breaks(self):
        """The 'breaks' command."""
        probes = cnvlib.read("formats/amplicon.cnr")
        segs = segmentation.do_segmentation("formats/amplicon.cnr", False,
                                            "haar")
        rows = commands.do_breaks(probes, segs, 4)
        self.assertTrue(len(rows) > 0)

    def test_call(self):
        """The 'call' command."""
        # Methods: clonal, threshold
        tr_cns = cnvlib.read("formats/tr95t.cns")
        tr_thresh = commands.do_call(tr_cns, "threshold",
                            is_reference_male=True, is_sample_female=True)
        self.assertEqual(len(tr_cns), len(tr_thresh))
        tr_clonal = commands.do_call(tr_cns, "clonal",
                            purity=.65,
                            is_reference_male=True, is_sample_female=True)
        self.assertEqual(len(tr_cns), len(tr_clonal))
        cl_cns = cnvlib.read("formats/cl_seq.cns")
        cl_thresh = commands.do_call(cl_cns, "threshold",
                            thresholds=np.log2((np.arange(12) + .5) / 6.),
                            is_reference_male=True, is_sample_female=True)
        self.assertEqual(len(cl_cns), len(cl_thresh))
        cl_clonal = commands.do_call(cl_cns, "clonal",
                            ploidy=6, purity=.99,
                            is_reference_male=True, is_sample_female=True)
        self.assertEqual(len(cl_cns), len(cl_clonal))

    def test_call_gender(self):
        """Test each 'call' method on allosomes."""
        for (fname, sample_is_f, ref_is_m, chr1_expect, chrx_expect, chry_expect
            ) in (
                ("formats/f-on-f.cns", True, False, 0, 0, 0),
                ("formats/f-on-m.cns", True, True, 0.585, 1, -9.97),
                ("formats/m-on-f.cns", False, False, 0, -1, 0),
                ("formats/m-on-m.cns", False, True, 0, 0, 0),
            ):
            cns = cnvlib.read(fname)
            chr1_idx = (cns.chromosome == 'chr1')
            chrx_idx = (cns.chromosome == 'chrX')
            chry_idx = (cns.chromosome == 'chrY')
            def test_chrom_means(segments):
                self.assertAlmostEqual(chr1_expect,
                                       segments['log2'][chr1_idx].mean(), 2)
                self.assertAlmostEqual(chrx_expect,
                                       segments['log2'][chrx_idx].mean(), 2)
                self.assertAlmostEqual(chry_expect,
                                       segments['log2'][chry_idx].mean(), 2)

            # Call threshold
            cns_thresh = commands.do_call(cns, "threshold",
                                 is_reference_male=ref_is_m,
                                 is_sample_female=sample_is_f)
            test_chrom_means(cns_thresh)
            # Call clonal pure
            cns_clone = commands.do_call(cns, "clonal",
                                is_reference_male=ref_is_m,
                                is_sample_female=sample_is_f)
            test_chrom_means(cns_clone)
            # Call clonal barely-mixed
            cns_p99 = commands.do_call(cns, "clonal", purity=0.99,
                              is_reference_male=ref_is_m,
                              is_sample_female=sample_is_f)
            test_chrom_means(cns_p99)

    def test_export(self):
        """Run the 'export' command with each format."""
        # SEG
        seg_rows = export.export_seg(["formats/tr95t.cns"])
        self.assertTrue(len(seg_rows) > 0)
        seg2_rows = export.export_seg(["formats/tr95t.cns",
                                       "formats/cl_seq.cns"])
        self.assertTrue(len(seg2_rows) > len(seg_rows))
        # THetA2
        _header, theta_rows = export.export_theta("formats/tr95t.cns",
                                                  "formats/reference-tr.cnn")
        self.assertTrue(len(theta_rows) > 0)
        # VCF
        tr_cns = cnvlib.read("formats/tr95t.cns")
        _header, tr_vcf_body = export.export_vcf(tr_cns, 2, True, True)
        self.assertTrue(0 < len(tr_vcf_body.splitlines()) < len(tr_cns))
        cl_cns = cnvlib.read("formats/cl_seq.cns")
        _header, cl_vcf_body = export.export_vcf(cl_cns, 6, True, True)
        self.assertTrue(0 < len(cl_vcf_body.splitlines()) < len(cl_cns))

    def test_gainloss(self):
        """The 'gainloss' command."""
        probes = cnvlib.read("formats/amplicon.cnr")
        rows = commands.do_gainloss(probes, male_reference=True)
        self.assertTrue(len(rows) > 0)
        segs = segmentation.do_segmentation("formats/amplicon.cnr", False,
                                            "haar")
        rows = commands.do_gainloss(probes, segs, True, 0.3, 4)
        self.assertTrue(len(rows) > 0)

    def test_metrics(self):
        """The 'metrics' command."""
        # TODO

    def test_reference(self):
        """The 'reference' command."""
        # Empty antitargets
        ref = commands.do_reference(["formats/amplicon.cnr"], ["formats/empty"])
        self.assertTrue(len(ref) > 0)
        # Empty antitargets, flat reference
        ref = commands.do_reference_flat("formats/amplicon.bed",
                                         "formats/empty")
        self.assertTrue(len(ref) > 0)

    def test_segment(self):
        """The 'segment' command."""
        # R methods are in another script
        cns = segmentation.do_segmentation("formats/amplicon.cnr", False,
                                           "haar")
        self.assertTrue(len(cns) > 0)

    def test_segmetrics(self):
        """The 'segmetrics' command."""
        # TODO

    def test_target(self):
        """The 'target' command."""
        # ENH: annotate w/ mini-refFlat
        commands.do_targets("formats/nv2_baits.interval_list", '/dev/null')
        commands.do_targets("formats/amplicon.bed", '/dev/null',
                            do_short_names=True, do_split=True, avg_size=100)



class OtherTests(unittest.TestCase):
    """Tests for other functionality."""

    def test_fix_edge(self):
        """Test the 'edge' bias correction calculations."""
        # With no gap, gain and loss should balance out
        # 1. Wide target, no secondary corrections triggered
        target_size = 600
        insert_size = 250
        loss = fix.edge_loss(target_size, insert_size)
        gain = fix.edge_gain(target_size, insert_size, 0)  # Adjacent
        gain *= 2  # Same on the other side
        self.assertAlmostEqual(loss, gain)
        # TODO - what tests make sense here?
        # 2. Trigger 'loss' correction (target_size < 2 * insert_size)
        # target_size = 600
        # self.assertAlmostEqual(fix.edge_loss(target_size, insert_size),
        #                        2 * fix.edge_gain(target_size, insert_size, 0))
        # 3. Trigger 'gain' correction (target_size + gap_size < insert_size)
        # target_size = 300
        # self.assertAlmostEqual(fix.edge_loss(target_size, insert_size),
        #                        2 * fix.edge_gain(target_size, insert_size, 0))

    # call
    # Test: convert_clonal(x, 1, 2) == convert_diploid(x)


# == helpers ==

def linecount(filename):
    i = 0
    with open(filename) as handle:
        for i, _line in enumerate(handle):
            pass
        return i + 1


if __name__ == '__main__':
    unittest.main(verbosity=2)
