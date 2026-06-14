#!/usr/bin/env python
"""Tests for the VariantArray class (cnvlib.vary)."""

import logging
import math
import unittest
import warnings

import pytest

logging.basicConfig(level=logging.ERROR, format="%(message)s")

import numpy as np
import pandas as pd

import cnvlib
from cnvlib import (
    cnary,
    fix,
    params,
    plots,
    reports,
    segfilters,
    segmentation,
    segmetrics,
    vary,
)
from skgenome import GenomicArray as GA
from skgenome import tabio


class VariantArrayTests(unittest.TestCase):
    """Tests for the VariantArray class."""

    def test_read(self):
        """Instantiate from a VCF file."""
        variants = tabio.read("formats/na12878_na12882_mix.vcf", "vcf")
        self.assertGreater(len(variants), 0)

    def test_baf_by_ranges_no_alt_freq_matches_ranges_length(self):
        """Without an alt_freq column, baf_by_ranges returns one NaN per *range*
        (not per variant), so callers assigning it as a column don't crash."""
        varr = vary.VariantArray(
            pd.DataFrame(
                {
                    "chromosome": ["chr1"] * 5,
                    "start": np.arange(5) * 100,
                    "end": np.arange(5) * 100 + 1,
                    "ref": ["A"] * 5,
                    "alt": ["G"] * 5,
                    "zygosity": np.full(5, 0.5),
                }
            )
        )
        ranges = cnary.CopyNumArray.from_rows(
            [["chr1", i * 100, i * 100 + 100, "-", 0.0] for i in range(3)],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )
        baf = varr.baf_by_ranges(ranges)
        self.assertEqual(len(baf), len(ranges))  # 3, not 5 (the variant count)
        self.assertTrue(baf.isna().all())

    def test_mirrored_baf_all_nan(self):
        """All-NaN frequencies mirror to all-NaN without warning (#407).

        np.nanmedian warns 'Mean of empty slice' on an all-NaN slice under
        older numpy; with pytest's filterwarnings=error that breaks the build
        (the 3.11-min job). _mirrored_baf must short-circuit before calling
        .median() on such a slice.
        """
        vals = pd.Series([np.nan, np.nan, np.nan])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            out = vary._mirrored_baf(vals)
        self.assertEqual(len(out), len(vals))
        self.assertTrue(out.isna().all())
        # A partially-NaN slice still mirrors normally (median over real values)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            mixed = vary._mirrored_baf(pd.Series([np.nan, 0.2, 0.3]))
        self.assertTrue(math.isnan(mixed.iloc[0]))
        # median 0.25 < 0.5 -> below-half branch (0.5 - shift); both already minor
        self.assertAlmostEqual(mixed.iloc[1], 0.2)
        self.assertAlmostEqual(mixed.iloc[2], 0.3)

    @staticmethod
    def _make_varr(**cols):
        """Build a VariantArray of `chr1` loci from extra-column keyword lists."""
        n = len(next(iter(cols.values())))
        base = {
            "chromosome": ["chr1"] * n,
            "start": list(range(100, 100 + n * 10, 10)),
            "end": list(range(101, 101 + n * 10, 10)),
            "ref": ["A"] * n,
            "alt": ["G"] * n,
        }
        base.update(cols)
        return vary.VariantArray(pd.DataFrame(base))

    @staticmethod
    def _one_bin(start, end):
        return cnary.CopyNumArray.from_rows(
            [["chr1", start, end, "-", 0.0]],
            columns=["chromosome", "start", "end", "gene", "log2"],
        )

    def test_tumor_boost_branches(self):
        # t<n -> 0.5*t/n ; t>=n -> 1 - 0.5*(1-t)/(1-n); equal -> stays
        out = vary._tumor_boost(np.array([0.3, 0.7, 0.5]), np.array([0.5, 0.5, 0.5]))
        np.testing.assert_allclose(out.to_numpy(), [0.3, 0.7, 0.5])

    def test_tumor_boost_requires_matched_normal(self):
        varr = self._make_varr(zygosity=[0.5, 0.5], alt_freq=[0.4, 0.6])
        with self.assertRaises(ValueError):
            varr.tumor_boost()

    def test_baf_counts_by_ranges_clamps_alt_count(self):
        # minor = min(alt, depth-alt); the 20 alt_count is clamped to depth (10)
        varr = self._make_varr(
            zygosity=[0.5, 0.5, 0.5],
            alt_count=[4.0, 6.0, 20.0],
            depth=[10.0, 10.0, 10.0],
        )
        minor, depths = varr.baf_counts_by_ranges(self._one_bin(90, 130))
        self.assertEqual(int(minor.iloc[0]), 8)  # 4 + 4 + 0
        self.assertEqual(int(depths.iloc[0]), 30)

    def test_baf_counts_by_ranges_missing_columns(self):
        varr = self._make_varr(zygosity=[0.5, 0.5])
        self.assertIsNone(varr.baf_counts_by_ranges(self._one_bin(90, 130)))

    def test_het_frac_by_ranges(self):
        # 2 of 4 SNVs heterozygous (zygosity not in {0.0, 1.0})
        varr = self._make_varr(zygosity=[0.5, 0.0, 0.5, 1.0])
        frac = varr.het_frac_by_ranges(self._one_bin(90, 200))
        self.assertAlmostEqual(frac.iloc[0], 0.5)

    def test_zygosity_from_freq(self):
        varr = self._make_varr(alt_freq=[0.0, 0.5, 1.0], zygosity=[0.5, 0.5, 0.5])
        out = varr.zygosity_from_freq(het_freq=0.25, hom_freq=0.95)
        self.assertEqual(list(out["zygosity"]), [0.0, 0.5, 1.0])
        # Original is unmodified (method copies)
        self.assertEqual(list(varr["zygosity"]), [0.5, 0.5, 0.5])


class ChrXHetDensityTests(unittest.TestCase):
    """Tests for ``vary.chrx_het_density_rejects_haploid`` (#341).

    A diploid X chromosome has SNP heterozygosity comparable to autosomes
    (~30% in panel-typical populations); a haploid X has essentially zero
    heterozygous calls modulo sequencing-error noise. The helper runs a
    one-sided binomial test of the observed chrX het count against a
    permissive haploid-X null (5% error-rate ceiling), rejecting when
    the observation cannot be reconciled with haploid X at the chosen
    alpha. This is the independent sex confirmer used in
    ``verify_sample_sex`` when a VCF is supplied.
    """

    def test_diploid_x_population_het_rate_rejects_haploid(self):
        """Realistic diploid-X panel (~30% het) confidently rejects haploid."""
        rejected, p = vary.chrx_het_density_rejects_haploid(
            n_chrx_total=50, n_chrx_het=15
        )
        self.assertTrue(rejected)
        self.assertIsNotNone(p)
        self.assertLess(p, 0.001)

    def test_haploid_x_with_sequencing_noise_does_not_reject(self):
        """True haploid X with a sprinkling of error-rate hets is not rejected."""
        # 50 total chrX SNPs, 1 het (matches the 5% ceiling, well within noise)
        rejected, _ = vary.chrx_het_density_rejects_haploid(
            n_chrx_total=50, n_chrx_het=1
        )
        self.assertFalse(rejected)

    def test_underpowered_below_min_snps(self):
        """Fewer than the min-snps floor returns (False, None) regardless."""
        # Even with 100% het rate, can't reject with only 5 SNPs.
        rejected, p = vary.chrx_het_density_rejects_haploid(
            n_chrx_total=5, n_chrx_het=5
        )
        self.assertFalse(rejected)
        self.assertIsNone(p)

    def test_borderline_n10_needs_strong_signal(self):
        """With only 10 SNPs the test demands a very high het rate to reject.

        Per calibration (binom.sf at alpha=0.001, p_null=0.05): n=10 needs
        k>=5 (50% het rate) to reject; k=4 (40%) is not enough.
        """
        rejected_at_5, _ = vary.chrx_het_density_rejects_haploid(
            n_chrx_total=10, n_chrx_het=5
        )
        rejected_at_4, _ = vary.chrx_het_density_rejects_haploid(
            n_chrx_total=10, n_chrx_het=4
        )
        self.assertTrue(rejected_at_5)
        self.assertFalse(rejected_at_4)

    def test_zero_hets_never_rejects(self):
        """No observed hets → cannot reject haploid (the null prediction)."""
        rejected, _ = vary.chrx_het_density_rejects_haploid(
            n_chrx_total=100, n_chrx_het=0
        )
        self.assertFalse(rejected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
