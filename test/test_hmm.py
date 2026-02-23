#!/usr/bin/env python
"""Unit tests for the HMM segmentation module."""

import logging
import unittest

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose, assert_array_equal
from scipy.special import betaln, gammaln
from scipy.stats import betabinom as sp_betabinom

from cnvlib.segmentation.hmm import (
    BETABINOM_RHO,
    _batched_emission_baf,
    _betabinom_logpmf,
    _expected_log2,
    _expected_minor_freq,
    enumerate_states,
    log_emission_probs,
    variants_in_segment,
)
from cnvlib.vary import VariantArray

logging.basicConfig(level=logging.ERROR, format="%(message)s")


class BetaBinomLogPMFTests(unittest.TestCase):
    """Tests for the beta-binomial log PMF helper."""

    def test_against_scipy(self):
        """Matches scipy.stats.betabinom for a range of inputs."""
        k = np.array([0, 5, 10, 25, 50])
        n = np.array([100, 100, 100, 100, 100])
        a = np.array([50.0, 50.0, 50.0, 50.0, 50.0])
        b = np.array([50.0, 50.0, 50.0, 50.0, 50.0])
        result = _betabinom_logpmf(k, n, a, b)
        expected = sp_betabinom.logpmf(k, n, a, b)
        assert_allclose(result, expected, rtol=1e-10)

    def test_zero_depth(self):
        """Returns 0.0 when n=0 (no observation)."""
        k = np.array([0.0])
        n = np.array([0.0])
        a = np.array([100.0])
        b = np.array([100.0])
        result = _betabinom_logpmf(k, n, a, b)
        assert_allclose(result, 0.0, atol=1e-15)

    def test_broadcasting(self):
        """Broadcasts correctly for (T,1) x (1,N) shapes."""
        k = np.array([10.0, 20.0, 30.0])[:, np.newaxis]  # (3, 1)
        n = np.array([100.0, 100.0, 100.0])[:, np.newaxis]  # (3, 1)
        a = np.array([50.0, 100.0])[np.newaxis, :]  # (1, 2)
        b = np.array([150.0, 100.0])[np.newaxis, :]  # (1, 2)
        result = _betabinom_logpmf(k, n, a, b)
        self.assertEqual(result.shape, (3, 2))
        # All values should be finite negative log-probs
        self.assertTrue(np.all(np.isfinite(result)))
        self.assertTrue(np.all(result <= 0))

    def test_extreme_alpha_beta(self):
        """Handles floor alpha/beta of 0.5 without NaN."""
        k = np.array([0.0, 1.0])
        n = np.array([10.0, 10.0])
        a = np.array([0.5, 0.5])
        b = np.array([200.0, 200.0])
        result = _betabinom_logpmf(k, n, a, b)
        self.assertTrue(np.all(np.isfinite(result)))

    def test_high_depth_concentrates(self):
        """Higher depth gives more peaked distribution (lower entropy)."""
        # With rho=200, alpha=100, beta=100 → expected freq = 0.5
        a, b = 100.0, 100.0
        # Low depth: k=5 of n=10
        lp_low = _betabinom_logpmf(
            np.array([5.0]),
            np.array([10.0]),
            np.array([a]),
            np.array([b]),
        ).item()
        # High depth: k=50 of n=100
        lp_high = _betabinom_logpmf(
            np.array([50.0]),
            np.array([100.0]),
            np.array([a]),
            np.array([b]),
        ).item()
        # The mode log-prob should be higher (less negative) for low depth
        # because the distribution is more spread out
        self.assertGreater(lp_low, lp_high)


class ExpectedMinorFreqTests(unittest.TestCase):
    """Tests for _expected_minor_freq."""

    def test_diploid_het(self):
        """Diploid heterozygous (CN=2, minor=1) at any purity → 0.5."""
        minor = np.array([1.0])
        cn = np.array([2.0])
        result = _expected_minor_freq(minor, cn, purity=0.8)
        assert_allclose(result, 0.5)

    def test_pure_loh(self):
        """CN=2, minor=0 at purity=1.0 → freq=0.0."""
        minor = np.array([0.0])
        cn = np.array([2.0])
        result = _expected_minor_freq(minor, cn, purity=1.0)
        assert_allclose(result, 0.0)

    def test_cn_zero(self):
        """CN=0 returns 0.5 (no allelic information)."""
        minor = np.array([0.0])
        cn = np.array([0.0])
        result = _expected_minor_freq(minor, cn, purity=1.0)
        assert_allclose(result, 0.5)

    def test_impure_loh(self):
        """LOH at 50% purity: (0.5*0 + 0.5*1) / (0.5*2 + 0.5*2) = 0.25."""
        minor = np.array([0.0])
        cn = np.array([2.0])
        result = _expected_minor_freq(minor, cn, purity=0.5)
        assert_allclose(result, 0.25)

    def test_values_at_most_half(self):
        """All expected minor freqs should be <= 0.5."""
        states = enumerate_states(6, include_baf=True)
        minor = np.array(
            [s[1] if s[1] is not None else 0 for s in states], dtype=np.float64
        )
        cn = np.array([s[0] for s in states], dtype=np.float64)
        for purity in [0.2, 0.5, 0.8, 1.0]:
            freqs = _expected_minor_freq(minor, cn, purity)
            self.assertTrue(
                np.all(freqs <= 0.5 + 1e-10), f"freq > 0.5 at purity={purity}: {freqs}"
            )


class LogEmissionProbsTests(unittest.TestCase):
    """Tests for log_emission_probs with beta-binomial BAF."""

    def setUp(self):
        self.states = enumerate_states(4, include_baf=True)
        self.n_states = len(self.states)

    def test_log2_only(self):
        """Works without BAF data (minor_counts=None, depths=None)."""
        log2_obs = np.array([0.0, -0.5, 0.3], dtype=np.float64)
        result = log_emission_probs(
            log2_obs,
            None,
            None,
            self.states,
            1.0,
            2.0,
            0.2,
        )
        self.assertEqual(result.shape, (3, self.n_states))
        self.assertTrue(np.all(np.isfinite(result)))

    def test_with_baf_counts(self):
        """BAF data adds to emission probabilities."""
        log2_obs = np.array([0.0, -0.5, 0.3], dtype=np.float64)
        minor_counts = np.array([50.0, 10.0, 30.0], dtype=np.float64)
        depths = np.array([100.0, 100.0, 100.0], dtype=np.float64)

        result_no_baf = log_emission_probs(
            log2_obs,
            None,
            None,
            self.states,
            0.8,
            2.0,
            0.2,
        )
        result_with_baf = log_emission_probs(
            log2_obs,
            minor_counts,
            depths,
            self.states,
            0.8,
            2.0,
            0.2,
        )
        self.assertEqual(result_with_baf.shape, result_no_baf.shape)
        # BAF should change the probabilities
        self.assertFalse(np.allclose(result_no_baf, result_with_baf))

    def test_zero_depth_bins_neutral(self):
        """Bins with zero depth contribute nothing to BAF emission."""
        log2_obs = np.array([0.0, 0.0], dtype=np.float64)
        minor_counts = np.array([0.0, 0.0], dtype=np.float64)
        depths = np.array([0.0, 0.0], dtype=np.float64)

        result_with_baf = log_emission_probs(
            log2_obs,
            minor_counts,
            depths,
            self.states,
            0.8,
            2.0,
            0.2,
        )
        result_no_baf = log_emission_probs(
            log2_obs,
            None,
            None,
            self.states,
            0.8,
            2.0,
            0.2,
        )
        assert_allclose(result_with_baf, result_no_baf)

    def test_het_baf_favors_diploid_het(self):
        """Balanced allele counts should favor CN=2, minor=1 state."""
        log2_obs = np.array([0.0], dtype=np.float64)
        minor_counts = np.array([50.0], dtype=np.float64)
        depths = np.array([100.0], dtype=np.float64)

        result = log_emission_probs(
            log2_obs,
            minor_counts,
            depths,
            self.states,
            1.0,
            2.0,
            0.2,
        )
        # Find the (CN=2, minor=1) state
        het_idx = next(i for i, s in enumerate(self.states) if s == (2, 1))
        # It should have the highest emission probability
        self.assertEqual(result[0].argmax(), het_idx)


class BatchedEmissionBAFTests(unittest.TestCase):
    """Tests for _batched_emission_baf."""

    def test_shape(self):
        """Output shape is (G, T, N) = (n_purities * n_ploidies, T, N)."""
        states = enumerate_states(4, include_baf=True)
        cn_arr = np.array([s[0] for s in states], dtype=np.float64)
        minor_arr = np.array(
            [s[1] if s[1] is not None else 0 for s in states], dtype=np.float64
        )
        purities = np.array([0.5, 0.8, 1.0])
        n_ploidies = 2
        T = 10
        mc = np.random.default_rng(42).integers(0, 50, size=T).astype(np.float64)
        dp = mc + np.random.default_rng(43).integers(0, 50, size=T).astype(np.float64)

        result = _batched_emission_baf(
            mc,
            dp,
            cn_arr,
            minor_arr,
            purities,
            n_ploidies,
            rho=BETABINOM_RHO,
        )
        G = len(purities) * n_ploidies
        N = len(states)
        self.assertEqual(result.shape, (G, T, N))
        self.assertTrue(np.all(np.isfinite(result)))

    def test_zero_depth_zeroed(self):
        """Zero-depth bins produce 0.0 emission (no BAF information)."""
        states = enumerate_states(4, include_baf=True)
        cn_arr = np.array([s[0] for s in states], dtype=np.float64)
        minor_arr = np.array(
            [s[1] if s[1] is not None else 0 for s in states], dtype=np.float64
        )
        purities = np.array([0.8])
        mc = np.array([0.0, 25.0, 0.0])
        dp = np.array([0.0, 100.0, 0.0])

        result = _batched_emission_baf(
            mc,
            dp,
            cn_arr,
            minor_arr,
            purities,
            1,
            rho=BETABINOM_RHO,
        )
        # Bins 0 and 2 (depth=0) should be all zeros
        assert_array_equal(result[0, 0, :], 0.0)
        assert_array_equal(result[0, 2, :], 0.0)
        # Bin 1 (depth=100) should be non-zero
        self.assertTrue(np.any(result[0, 1, :] != 0.0))


def _make_variant_array(chroms, starts, ends, alt_freqs):
    """Create a VariantArray with the given positions and allele frequencies."""
    n = len(starts)
    return VariantArray(
        pd.DataFrame(
            {
                "chromosome": chroms,
                "start": starts,
                "end": ends,
                "ref": ["C"] * n,
                "alt": ["T"] * n,
                "alt_freq": alt_freqs,
            }
        )
    )


def _make_segment(chrom, start, end, log2=0.0, gene="-", probes=100):
    """Create a namedtuple-like segment row."""
    return pd.Series(
        {
            "chromosome": chrom,
            "start": start,
            "end": end,
            "log2": log2,
            "gene": gene,
            "probes": probes,
        }
    )


class VariantsInSegmentTests(unittest.TestCase):
    """Tests for variants_in_segment (BAF re-segmentation for non-HMM methods)."""

    def test_too_few_variants_returns_original_segment(self):
        """With fewer variants than MIN_VARIANTS_THRESHOLD, returns the
        original segment unchanged."""
        varr = _make_variant_array(
            ["chr1"] * 10,
            list(range(100, 200, 10)),
            list(range(101, 201, 10)),
            [0.5] * 10,
        )
        segment = _make_segment("chr1", 0, 1000, log2=-0.3)
        result = variants_in_segment(varr, segment)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["chromosome"], "chr1")
        self.assertEqual(result.iloc[0]["start"], 0)
        self.assertEqual(result.iloc[0]["end"], 1000)
        self.assertAlmostEqual(result.iloc[0]["log2"], -0.3)

    def test_uniform_baf_returns_single_segment(self):
        """Uniform BAF near 0.5 should produce a single segment (no extra
        breakpoints)."""
        n = 200
        rng = np.random.default_rng(42)
        alt_freqs = rng.normal(0.5, 0.02, size=n).clip(0.01, 0.99)
        starts = np.arange(1000, 1000 + n * 100, 100)
        ends = starts + 1
        varr = _make_variant_array(
            ["chr1"] * n, starts.tolist(), ends.tolist(), alt_freqs.tolist()
        )
        segment = _make_segment("chr1", 0, 1000 + n * 100)
        result = variants_in_segment(varr, segment, min_variants=50)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["start"], 0)
        self.assertEqual(result.iloc[0]["end"], 1000 + n * 100)

    def test_shifted_baf_creates_breakpoint(self):
        """A region with BAF ~0.5 followed by BAF ~0.67 should produce at
        least 2 segments."""
        rng = np.random.default_rng(123)
        n_half = 100
        baf_neutral = rng.normal(0.5, 0.02, size=n_half).clip(0.01, 0.99)
        baf_shifted = rng.normal(0.67, 0.02, size=n_half).clip(0.01, 0.99)
        alt_freqs = np.concatenate([baf_neutral, baf_shifted])
        n = len(alt_freqs)
        starts = np.arange(1000, 1000 + n * 100, 100)
        ends = starts + 1
        varr = _make_variant_array(
            ["chr1"] * n, starts.tolist(), ends.tolist(), alt_freqs.tolist()
        )
        segment = _make_segment("chr1", 0, 1000 + n * 100, probes=n)
        result = variants_in_segment(varr, segment, min_variants=50)
        self.assertGreaterEqual(len(result), 2)
        # All segments should have valid coordinates
        self.assertTrue((result["start"] < result["end"]).all())
        # First segment starts at original start, last ends at original end
        self.assertEqual(result.iloc[0]["start"], 0)
        self.assertEqual(result.iloc[-1]["end"], 1000 + n * 100)

    def test_output_columns(self):
        """Output DataFrame has the expected columns."""
        varr = _make_variant_array(
            ["chr1"] * 5,
            list(range(100, 600, 100)),
            list(range(101, 601, 100)),
            [0.5] * 5,
        )
        segment = _make_segment("chr1", 0, 1000)
        result = variants_in_segment(varr, segment)
        for col in ("chromosome", "start", "end", "gene", "log2", "probes"):
            self.assertIn(col, result.columns)

    def test_preserves_segment_metadata(self):
        """Gene name and log2 from the input segment are preserved."""
        varr = _make_variant_array(
            ["chr1"] * 5,
            list(range(100, 600, 100)),
            list(range(101, 601, 100)),
            [0.5] * 5,
        )
        segment = _make_segment("chr1", 0, 1000, log2=-0.5, gene="TP53")
        result = variants_in_segment(varr, segment)
        self.assertEqual(result.iloc[0]["gene"], "TP53")
        self.assertAlmostEqual(result.iloc[0]["log2"], -0.5)


if __name__ == "__main__":
    unittest.main()
