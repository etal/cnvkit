#!/usr/bin/env python
"""Unit tests for the 'genome' sub-package."""
from __future__ import absolute_import, division, print_function

import unittest

from cnvlib import tabio


class GaryTests(unittest.TestCase):

    def setUp(self):
        self.ex_cnr = tabio.read_cna('formats/reference-tr.cnn')

    def test_empty(self):
        """Instantiate from an empty file."""
        garr = tabio.read("formats/empty")
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

    def test_autosomes(self):
        """Test selection of autosomes."""
        len_all = len(self.ex_cnr)
        len_x = (self.ex_cnr.chromosome == 'chrX').sum()
        len_y = (self.ex_cnr.chromosome == 'chrY').sum()
        auto = self.ex_cnr.autosomes()
        self.assertEqual(len(auto), len_all - len_x - len_y)
        autox = self.ex_cnr.autosomes(also='chrX')
        self.assertEqual(len(autox), len_all - len_y)
        autoy = self.ex_cnr.autosomes(also=['chrY'])
        self.assertEqual(len(autoy), len_all - len_x)
        autoxy = self.ex_cnr.autosomes(also=['chrX', 'chrY'])
        self.assertEqual(len(autoxy), len_all)

    def test_by_chromosome(self):
        for fname in ("formats/amplicon.cnr", "formats/cl_seq.cns"):
            cnarr = tabio.read_cna(fname)
            row_count = 0
            for _chrom, rows in cnarr.by_chromosome():
                row_count += len(rows)
            self.assertEqual(row_count, len(cnarr))

    # def test_concat(self):

    def test_filter(self):
        """Test sugary selection of a subset of the data array."""
        num_bg_rows = len(self.ex_cnr[self.ex_cnr['gene'] == 'Background'])
        self.assertEqual(len(self.ex_cnr.filter(gene='Background')),
                         num_bg_rows)
        selector = lambda row: row['gene'] == 'Background'
        self.assertEqual(len(self.ex_cnr.filter(selector)), num_bg_rows)

    def test_ranges(self):
        """Test range methods: by_ranges, in_range, in_ranges."""
        cnarr = tabio.read_cna("formats/amplicon.cnr")
        segarr = tabio.read_cna("formats/amplicon.cns")
        chrom_segarr = dict(segarr.by_chromosome())
        for chrom, subarr in cnarr.by_chromosome():
            count_segs = 0
            count_bins = 0
            subsegarr = chrom_segarr[chrom]
            for count_segs, (seg, bins) in enumerate(subarr.by_ranges(subsegarr)):
                count_bins += len(bins)
                self.assertEqual(seg.probes, len(bins))
                self.assertEqual(len(bins), len(
                    cnarr.in_range(seg.chromosome, seg.start, seg.end,
                                   mode='outer')))
                self.assertEqual(len(bins), len(
                    cnarr.in_range(seg.chromosome, seg.start, seg.end,
                                   mode='trim')))
            self.assertEqual(len(subsegarr), count_segs + 1)
            self.assertEqual(len(subarr), count_bins)
            self.assertEqual(len(subarr), len(
                cnarr.in_ranges(chrom, subsegarr['start'], subsegarr['end'],
                                mode="outer")))
            self.assertEqual(len(subarr), len(
                subarr.in_ranges(starts=subsegarr['start'],
                                 ends=subsegarr['end'], mode="outer")))
            self.assertEqual(len(subarr), len(
                cnarr.in_ranges(chrom, subsegarr['start'], subsegarr['end'],
                                mode="trim")))
            self.assertEqual(len(subarr), len(
                subarr.in_ranges(starts=subsegarr['start'],
                                 ends=subsegarr['end'], mode="trim")))

    def test_resize_ranges(self):
        """Test resizing bins."""
        baits_fname = 'formats/nv2_baits.interval_list'
        chrom_sizes = {'chr1': 249250621,
                       'chr2': 243199373,
                       'chr3': 198022430,
                       'chr4': 191154276,
                       'chr5': 180915260,
                       'chr6': 171115067,
                       'chr7': 159138663,
                       'chr8': 146364022,
                       'chr9': 141213431,
                       'chr10': 135534747,
                       'chr11': 135006516,
                       'chr12': 133851895,
                       'chr13': 115169878,
                       'chr14': 107349540,
                       'chr15': 102531392,
                       'chr16': 90354753,
                       'chr17': 81195210,
                       'chr18': 78077248,
                       'chr19': 59128983,
                       'chr20': 63025520,
                       'chr21': 48129895,
                       'chr22': 51304566,
                       'chrX': 155270560,
                       'chrY': 59373566}
        bins = tabio.read(baits_fname, 'interval')
        for chrom, arr in bins.resize_ranges(1e7, chrom_sizes).by_chromosome():
            self.assertLessEqual(0, arr.start.min())
            self.assertLessEqual(arr.end.max(), chrom_sizes[chrom])

    def test_shuffle_sort(self):
        """Test shuffling and re-sorting the data array."""
        orig_cvg = tuple(self.ex_cnr['log2'][:10])
        self.assertEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)
        self.ex_cnr.shuffle()
        self.assertNotEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)
        self.ex_cnr.sort()
        self.assertEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)



if __name__ == '__main__':
    unittest.main(verbosity=2)
