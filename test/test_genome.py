#!/usr/bin/env python
"""Unit tests for the 'genome' sub-package."""
from __future__ import absolute_import, division, print_function
import unittest

import pandas as pd

from cnvlib import tabio
from cnvlib.genome import GenomicArray as GA


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

    def test_total_range_size(self):
        """Test total region coverage calculation."""
        for fname, area in (
            ('formats/empty', 0),
            ('formats/duke-my.bed', 103),
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

    # Likely overlapping gene models
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
                         #  '\n'.join(["Got:", str(result.data),
                         #             "Expected:", str(expect.data)])
                        )
        for col in expect.data.columns:
            self.assertTrue((expect[col] == result[col]).all(),
                            "Col '{}' differs:\nExpect:\n{}\nGot:\n{}"
                            #  .format(col, expect[col], result[col]))
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
        # =A=========================
        #   =B=======   =D===   =E======
        #      =C=
        # 1 3  5 6  8   11 15   19 20 23 <- coordinates

        selections = self._from_intervals([
            (1, 8, ''),
            (4, 10, ''),
            (8, 19, ''),
            (11, 20, ''),
            (21, 22, ''),
        ])

        expectations = {
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

        for mode in ('outer', 'trim', 'inner'):
            grouped_results = self.regions_1.by_ranges(selections, mode=mode)
            for (_coord, result), expect in zip(grouped_results,
                                                expectations[mode]):
                self._compare_regions(result, self._from_intervals(expect))

        # TODO region_coords_2
        # =A=============================
        #   =B==  =C==     =E==  =G==
        #               =D=========================
        #                  =F==  =H==       =I==
        # 3 5  8  11 14 17 19 22 25 28  32  36 39 42

        selections = self._from_intervals([
            (0, 1, ''),
        ])
        expectations = {
            'outer': (
                [],
            ),
            'trim': (
                [],
            ),

            'inner': (
                [],
            ),
        }
        for mode in ('outer', 'trim', 'inner'):
            grouped_results = self.regions_2.by_ranges(selections, mode=mode)
            for (_coord, result), expect in zip(grouped_results,
                                                expectations[mode]):
                self._compare_regions(result, self._from_intervals(expect))



if __name__ == '__main__':
    unittest.main(verbosity=2)
