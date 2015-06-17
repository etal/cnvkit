#!/usr/bin/env python

"""Unit tests for the CNVkit library, cnvlib."""

import unittest

import numpy
# from Bio._py3k import StringIO

import cnvlib
# Import all modules as a smoke test
from cnvlib import (antitarget, commands, core, coverage, diagram, export, fix,
                    importers, metrics, ngfrills, params, plots, cnarray,
                    reference, reports, segmentation, smoothing)

EX_CNR = 'reference-tr.cnn'


class CNATests(unittest.TestCase):
    """Tests for the CopyNumArray class."""

    def setUp(self):
        self.ex_cnr = cnvlib.read(EX_CNR)

    def test_basic(self):
        """Test basic container functionality and magic methods."""
        # Length
        self.assertEqual(len(self.ex_cnr), 27526)
        # Equality
        same = cnvlib.read(EX_CNR)
        self.assertEqual(self.ex_cnr, same)
        # Item access
        orig = self.ex_cnr[0]
        self.ex_cnr[0] = orig
        self.ex_cnr[3:4] = self.ex_cnr[3:4]
        self.ex_cnr[6:10] = self.ex_cnr[6:10]
        self.assertEqual(tuple(self.ex_cnr[0]), tuple(same[0]))
        self.assertEqual(self.ex_cnr[3:6], same[3:6])

    def test_iter(self):
        """Test iteration."""
        rows = iter(self.ex_cnr)
        firstrow = next(rows)
        self.assertEqual(tuple(firstrow), tuple(self.ex_cnr[0]))
        i = 0
        for i, row in enumerate(rows):
            pass
        self.assertEqual(i + 2, len(self.ex_cnr))

    # def test_by_bin(self):
    # def test_by_chromosome(self):
    # def test_by_gene(self):
    # def test_by_segment(self):

    def test_center_all(self):
        """Test median-recentering."""
        chr1 = self.ex_cnr.in_range('chr1')
        self.assertAlmostEqual(0, numpy.median(chr1['log2']), places=1)
        chr1.center_all()
        orig_chr1_cvg = numpy.median(chr1['log2'])
        self.assertAlmostEqual(0, orig_chr1_cvg)
        chr1plus2 = chr1.copy()
        chr1plus2['log2'] += 2.0
        chr1plus2.center_all()
        self.assertAlmostEqual(numpy.median(chr1plus2['log2']), orig_chr1_cvg)

    def test_copy(self):
        """Test creation of an independent copy of the object."""
        dupe = self.ex_cnr.copy()
        self.assertEqual(tuple(self.ex_cnr[3]), tuple(dupe[3]))
        self.ex_cnr[3, 'log2'] = -10.0
        self.assertNotEqual(tuple(self.ex_cnr[3]), tuple(dupe[3]))

    def test_drop_extra_columns(self):
        """Test removal of optional 'gc' column."""
        self.assertTrue('gc' in self.ex_cnr)
        cleaned = self.ex_cnr.drop_extra_columns()
        self.assertTrue('gc' not in cleaned)
        self.assertTrue((cleaned['log2'] == self.ex_cnr['log2']).all())

    # def test_extend(self):
    # def test_in_range(self):

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

    # def test_squash_genes(self):


class ImporterTests(unittest.TestCase):
    """Tests for importers functionality."""

    def test_import_picard(self):
        """Test loading a Picard targetcoverage file."""
        fname = 'picard/p2-5_5.antitargetcoverage.csv'
        cna = importers.load_targetcoverage_csv(fname)
        self.assertTrue(len(cna) > 1)


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



if __name__ == '__main__':
    unittest.main(verbosity=2)
