#!/usr/bin/env python
"""Unit tests for the 'genome' sub-package."""
from __future__ import absolute_import, division, print_function
import random
import unittest

import numpy as np
import pandas as pd

from cnvlib import read
from skgenome import chromsort, rangelabel
from skgenome import tabio, GenomicArray as GA


class GaryTests(unittest.TestCase):

    def setUp(self):
        self.ex_cnr = read('formats/reference-tr.cnn')

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
            cnarr = read(fname)
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

    def test_ranges_by_in(self):
        """Test range methods: by_ranges, in_range, in_ranges."""
        cnarr = read("formats/amplicon.cnr")
        segarr = read("formats/amplicon.cns")
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

    def test_ranges_into(self):
        cnarr = read("formats/amplicon.cnr")
        segarr = read("formats/amplicon.cns")
        seg_genes = cnarr.into_ranges(segarr, 'gene', '-')
        self.assertEqual(len(seg_genes), len(segarr))
        # With a VCF
        varr = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        seg_baf = varr.into_ranges(segarr, 'alt_freq', np.nan, np.nanmedian)
        self.assertEqual(len(seg_baf), len(segarr))
        cna_baf = varr.into_ranges(cnarr, 'alt_freq', 0.0, np.max)
        self.assertEqual(len(cna_baf), len(cnarr))
        # Edge cases
        mtarr = tabio.read("formats/empty")
        segarr.into_ranges(mtarr, 'start', 0, int)
        mtarr.into_ranges(segarr, 'end', 0, 0)

    def test_ranges_resize(self):
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

    def test_total_range_size(self):
        """Test total region coverage calculation."""
        for fname, area in (
            ('formats/empty', 0),
            ('formats/my-targets.bed', 103),
            ('formats/dac-my.bed', 148),
            ('formats/example.gff', 7951),
            ('formats/refflat-mini.txt', 719715),
        ):
            regions = tabio.read_auto(fname)
            self.assertEqual(regions.total_range_size(), area)

    def test_shuffle_sort(self):
        """Test shuffling and re-sorting the data array."""
        orig_cvg = tuple(self.ex_cnr['log2'][:10])
        self.assertEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)
        self.ex_cnr.shuffle()
        self.assertNotEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)
        self.ex_cnr.sort()
        self.assertEqual(tuple(self.ex_cnr['log2'][:10]), orig_cvg)



class IntervalTests(unittest.TestCase):
    """Interval arithmetic tests."""
    combiner = {'gene': lambda a: ''.join(a)}

    # Simple: nested, overlapping, & non-overlapping intervals
    # =A=========================
    #   =B=======   =D===   =E======
    #      =C=
    # 1 3  5 6  8   11 15   19 20 23 <- coordinates
    region_coords_1 = (
        (1,  20, 'A'),
        (3,  8,  'B'),
        (5,  6,  'C'),
        (11, 15, 'D'),
        (19, 23, 'E'),
    )

    # Semi-realistic: overlapping gene models
    # =A=============================
    #   =B==  =C==     =E==  =G==
    #               =D=========================
    #                  =F==  =H==       =I==
    # 3 5  8  11 14 17 19 22 25 28  32  36 39 42
    region_coords_2 = (
        (3,  32, 'A'),
        (5,   8, 'B'),
        (11, 14, 'C'),
        (17, 42, 'D'),
        (19, 22, 'E'),
        (19, 22, 'F'),
        (25, 28, 'G'),
        (25, 28, 'H'),
        (36, 39, 'I'),
    )

    @staticmethod
    def _from_intervals(coords):
        garr = GA(pd.DataFrame(list(coords),
                               columns=['start', 'end', 'gene'])
                  .assign(chromosome='chr0'))
        garr.sort_columns()
        return garr

    def _compare_regions(self, result, expect):
        self.assertEqual(expect.data.shape, result.data.shape,
                         '\n'.join(["Got:", str(result.data),
                                    "Expected:", str(expect.data)]))
        for col in expect.data.columns:
            self.assertTrue((expect[col] == result[col]).all(),
                            "Col '{}' differs:\nExpect:\n{}\nGot:\n{}"
                            .format(col, expect.data, result.data))

    def setUp(self):
        self.regions_1 = self._from_intervals(self.region_coords_1)
        self.regions_2 = self._from_intervals(self.region_coords_2)

    def test_flatten(self):
        flat_coords_1 = [
            (1,  3,  'A'),
            (3,  5,  'AB'),
            (5,  6,  'ABC'),
            (6,  8,  'AB'),
            (8,  11, 'A'),
            (11, 15, 'AD'),
            (15, 19, 'A'),
            (19, 20, 'AE'),
            (20, 23, 'E'),
        ]
        flat_coords_2 = [
            (3,  5,  'A'),
            (5,  8,  'AB'),
            (8,  11, 'A'),
            (11, 14, 'AC'),
            (14, 17, 'A'),
            (17, 19, 'AD'),
            (19, 22, 'ADEF'),
            (22, 25, 'AD'),
            (25, 28, 'ADGH'),
            (28, 32, 'AD'),
            (32, 36, 'D'),
            (36, 39, 'DI'),
            (39, 42, 'D'),
        ]
        for regions, flat_coords in [
            (self.regions_1, flat_coords_1),
            (self.regions_2, flat_coords_2),
        ]:
            result = regions.flatten(combine=self.combiner)
            expect = self._from_intervals(flat_coords)
            self._compare_regions(result, expect)

    def test_merge(self):
        merged_coords_1 = [(1, 23, 'ABCDE')]
        merged_coords_2 = [(3, 42, 'ABCDEFGHI')]
        for regions, merged_coords in [
            (self.regions_1, merged_coords_1),
            (self.regions_2, merged_coords_2),
        ]:
            result = regions.merge(combine=self.combiner)
            expect = self._from_intervals(merged_coords)
            self._compare_regions(result, expect)

    def test_intersect(self):
        selections1 = self._from_intervals([
            (1, 8, ''),
            (4, 10, ''),
            (8, 19, ''),
            (11, 20, ''),
            (21, 22, ''),
        ])
        expectations1 = {
            'outer': (
                # 1-8
                [(1,  20, 'A'),
                 (3,  8,  'B'),
                 (5,  6,  'C'),
                ],
                # 4-10
                [(1,  20, 'A'),
                 (3,  8,  'B'),
                 (5,  6,  'C'),
                ],
                # 8-19
                [(1,  20, 'A'),
                 (11, 15, 'D')],
                # 11-20
                [(1,  20, 'A'),
                 (11, 15, 'D'),
                 (19, 23, 'E')],
                # 21-22
                [(19, 23, 'E')],
            ),
            'trim': (
                # 1-8
                [(1,  8, 'A'),
                 (3,  8,  'B'),
                 (5,  6,  'C')],
                # 4-10
                [(4,  10, 'A'),
                 (4,  8,  'B'),
                 (5,  6,  'C')],
                # 8-19
                [(8,  19, 'A'),
                 (11, 15, 'D')],
                # 11-20
                [(11, 20, 'A'),
                 (11, 15, 'D'),
                 (19, 20, 'E')],
                # 21-22
                [(21, 22, 'E')],
            ),
            'inner': (
                # 1-8
                [(3,  8,  'B'),
                 (5,  6,  'C')],
                # 4-10
                [(5,  6,  'C')],
                # 8-19
                [(11, 15, 'D')],
                # 11-20
                [(11, 15, 'D')],
                # 21-22
                [],
            ),
        }

        selections2 = self._from_intervals([
            (0, 1, ''),
            (5, 14, ''),
            (16, 45, ''),
            (18, 37, ''),
            (19, 25, ''),
            (29, 31, ''),
            (34, 39, ''),
        ])
        expectations2 = {
            'outer': (
                # 0-1
                [],
                # 5-14
                [(3,  32, 'A'),
                 (5,   8, 'B'),
                 (11, 14, 'C')],
                # 16-45
                [(3,  32, 'A'),
                 (17, 42, 'D'),
                 (19, 22, 'E'),
                 (19, 22, 'F'),
                 (25, 28, 'G'),
                 (25, 28, 'H'),
                 (36, 39, 'I')],
                # 18-37
                [(3,  32, 'A'),
                 (17, 42, 'D'),
                 (19, 22, 'E'),
                 (19, 22, 'F'),
                 (25, 28, 'G'),
                 (25, 28, 'H'),
                 (36, 39, 'I')],
                # 19-25
                [(3,  32, 'A'),
                 (17, 42, 'D'),
                 (19, 22, 'E'),
                 (19, 22, 'F')],
                # 29-31
                [(3,  32, 'A'),
                 (17, 42, 'D')],
                # 34-39
                [(17, 42, 'D'),
                 (36, 39, 'I')],
            ),
            'trim': (
                # 0-1
                [],
                # 5-14
                [(5,  14, 'A'),
                 (5,   8, 'B'),
                 (11, 14, 'C')],
                # 16-45
                [(16, 32, 'A'),
                 (17, 42, 'D'),
                 (19, 22, 'E'),
                 (19, 22, 'F'),
                 (25, 28, 'G'),
                 (25, 28, 'H'),
                 (36, 39, 'I')],
                # 18-37
                [(18, 32, 'A'),
                 (18, 37, 'D'),
                 (19, 22, 'E'),
                 (19, 22, 'F'),
                 (25, 28, 'G'),
                 (25, 28, 'H'),
                 (36, 37, 'I')],
                # 19-25
                [(19, 25, 'A'),
                 (19, 25, 'D'),
                 (19, 22, 'E'),
                 (19, 22, 'F')],
                # 29-31
                [(29, 31, 'A'),
                 (29, 31, 'D')],
                # 34-39
                [(34, 39, 'D'),
                 (36, 39, 'I')],
            ),
            'inner': (
                # 0-1
                [],
                # 5-14
                [(5,   8, 'B'),
                 (11, 14, 'C')],
                # 16-45
                [(17, 42, 'D'),
                 (19, 22, 'E'),
                 (19, 22, 'F'),
                 (25, 28, 'G'),
                 (25, 28, 'H'),
                 (36, 39, 'I')],
                # 18-37
                [(19, 22, 'E'),
                 (19, 22, 'F'),
                 (25, 28, 'G'),
                 (25, 28, 'H')],
                # 19-25
                [(19, 22, 'E'),
                 (19, 22, 'F')],
                # 29-31
                [],
                # 34-39
                [(36, 39, 'I')],
            ),
        }

        for regions, selections, expectations in (
            (self.regions_1, selections1, expectations1),
            (self.regions_2, selections2, expectations2),
        ):
            for mode in ('outer', 'trim', 'inner'):
                grouped_results = regions.by_ranges(selections, mode=mode)
                for (_coord, result), expect in zip(grouped_results,
                                                    expectations[mode]):
                    self._compare_regions(result, self._from_intervals(expect))

    def test_subtract(self):
        # Test cases:
        #  | access: ====   ====   ====    ==========
        #  | target: ====  ====== ===   = ==  ==   ===
        #  | expect:                 ==     ==  ===
        #  | invert:       =    = =     = =          =
        #  |         1   5 78   2345 7 901234 6 8  1 34
        #  |        0         1         2         3

        access = self._from_intervals([
            (1, 5, ''),
            (8, 12, ''),
            (15, 19, ''),
            (23, 33, ''),
        ])
        target = self._from_intervals([
            (1, 5, ''),
            (7, 13, ''),
            (14, 17, ''),
            (20, 21, ''),
            (22, 24, ''),
            (26, 28, ''),
            (31, 34, ''),
        ])
        expect = self._from_intervals([
            (17, 19, ''),
            (24, 26, ''),
            (28, 31, ''),
        ])
        invert = self._from_intervals([
            (7, 8, ''),
            (12, 13, ''),
            (14, 15, ''),
            (20, 21, ''),
            (22, 23, ''),
            (33, 34, ''),
        ])

        result = access.subtract(target)
        self._compare_regions(result, expect)
        iresult = target.subtract(access)
        self._compare_regions(iresult, invert)



class OtherTests(unittest.TestCase):
    """Other small modules & functions in this sub-package."""

    def test_chromsort(self):
        labels_hg = ["chr1", "chr2", "chr10", "chr19", "chr20",
                     "chrX", "chrY", "chrM"]
        labels_grc = ["1", "2", "10", "19", "X", "Y", "MT"]
        for labels in (labels_hg, labels_grc):
            shuf = labels[:]
            random.shuffle(shuf)
            resorted = sorted(labels, key=chromsort.sorter_chrom)
            self.assertEqual(resorted, labels)

    def test_detect_big_chroms(self):
        sizes = [1, 20, 30]
        n_big, thresh = chromsort.detect_big_chroms(sizes)
        self.assertEqual(n_big, 2)
        self.assertEqual(thresh, 20)

    def test_rangelabel(self):
        row = rangelabel.from_label("chr1:123-456", keep_gene=False)
        self.assertEqual(tuple(row), ("chr1", 122, 456))
        label = rangelabel.to_label(row)
        self.assertEqual(label, "chr1:123-456")
        # unpack_range
        for region, expect in (
            ["chr1",                    ("chr1", None, None)],
            ["chr1:100-123",            ("chr1", 99, 123)],
            [("chr1", 100, 123),        ("chr1", 100, 123)],
            [("chr1", 100, 123, "A"),   ("chr1", 100, 123)],
        ):
            result = rangelabel.unpack_range(region)
            self.assertEqual(result, expect)



if __name__ == '__main__':
    unittest.main(verbosity=2)
