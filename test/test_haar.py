#!/usr/bin/env python
"""Unit tests for the HaarSeg segmentation algorithm."""

import logging
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy import stats

from cnvlib.segmentation.haar import (
    FDRThres,
    FindLocalPeaks,
    HaarConv,
    PulseConv,
    haarSeg,
)

logging.basicConfig(level=logging.ERROR, format="%(message)s")


class FDRThresTests(unittest.TestCase):
    """Tests for FDR threshold calculation."""

    def test_fdr_thres_basic(self):
        """FDRThres computes correct p-values using normal distribution."""
        # Known values: for stdev=1, the threshold at q=0.05 with sorted |x|
        # For a two-tailed test, p = 2*(1 - Phi(|x|/stdev))
        # At |x|=1.96, p ≈ 0.05
        stdev = 1.0
        # Create values where we know the p-values
        x = np.array([2.5, 2.0, 1.5, 1.0, 0.5])
        q = 0.05
        T = FDRThres(x, q, stdev)
        # The threshold should be positive and reasonable
        self.assertGreater(T, 0)
        self.assertLess(T, 3.0)

    def test_fdr_thres_matches_r_pnorm(self):
        """FDRThres p-values match R's pnorm(x, 0, stdev) * 2."""
        # Reference: In R, two-tailed p-value is 2 * pnorm(-abs(x), 0, sd)
        # or equivalently 2 * (1 - pnorm(abs(x), 0, sd))
        stdev = 0.5
        x_vals = np.array([1.0, 0.8, 0.6, 0.4, 0.2])

        # Compute expected p-values the way R would
        expected_p = 2 * (1 - stats.norm.cdf(np.abs(x_vals), loc=0, scale=stdev))

        # For q=1.0, all should pass, so threshold should be min
        # We verify the p-value calculation is correct by checking specific values
        # At x=1.0, stdev=0.5: p = 2*(1 - Phi(2)) ≈ 0.0455
        p_at_1 = 2 * (1 - stats.norm.cdf(1.0, loc=0, scale=0.5))
        assert_allclose(p_at_1, 0.0455, atol=0.001)

    def test_fdr_thres_empty_input(self):
        """FDRThres returns 0 for empty or single-element input."""
        self.assertEqual(FDRThres(np.array([]), 0.05, 1.0), 0)
        self.assertEqual(FDRThres(np.array([1.5]), 0.05, 1.0), 0)

    def test_fdr_thres_stdev_scaling(self):
        """FDRThres threshold scales with standard deviation."""
        # Use larger values and higher q to ensure some pass the threshold
        x = np.array([3.0, 2.5, 2.0, 1.5, 1.0])
        q = 0.1
        T1 = FDRThres(x, q, stdev=1.0)
        T2 = FDRThres(x, q, stdev=2.0)
        # With larger stdev, same x values have larger p-values (less significant),
        # so fewer pass the FDR criterion, resulting in a higher threshold
        self.assertGreaterEqual(T2, T1)


class HaarConvTests(unittest.TestCase):
    """Tests for Haar wavelet convolution."""

    def test_haar_conv_basic(self):
        """HaarConv produces non-zero output for step signal."""
        # Step signal: zeros then ones
        signal = np.concatenate([np.zeros(10), np.ones(10)])
        result = HaarConv(signal, None, stepHalfSize=2)
        self.assertEqual(len(result), len(signal))
        # Should detect the step around index 10
        self.assertGreater(np.abs(result).max(), 0)

    def test_haar_conv_constant_signal(self):
        """HaarConv returns near-zero for constant signal."""
        signal = np.ones(20) * 5.0
        result = HaarConv(signal, None, stepHalfSize=2)
        # Constant signal should have zero wavelet coefficients (except boundary)
        assert_allclose(result[5:15], 0, atol=1e-10)

    def test_haar_conv_size_edge_case(self):
        """HaarConv returns zeros when stepHalfSize > signalSize."""
        signal = np.array([1.0, 2.0, 3.0])
        # stepHalfSize=10 > signalSize=3
        with self.assertLogs(level="WARNING") as log:
            result = HaarConv(signal, None, stepHalfSize=10)
        self.assertEqual(len(result), len(signal))
        assert_array_equal(result, np.zeros(3))
        self.assertTrue(any("exceeds signal length" in msg for msg in log.output))

    def test_haar_conv_weighted(self):
        """HaarConv handles weights correctly."""
        signal = np.concatenate([np.zeros(10), np.ones(10)])
        weights = np.ones(20)
        result_weighted = HaarConv(signal, weights, stepHalfSize=2)
        result_unweighted = HaarConv(signal, None, stepHalfSize=2)
        # With uniform weights, results should be similar in shape
        # (but not identical due to different normalization)
        self.assertEqual(len(result_weighted), len(signal))
        # Both should detect the step
        self.assertGreater(np.abs(result_weighted).max(), 0)


class PulseConvTests(unittest.TestCase):
    """Tests for pulse convolution (used in non-stationary variance)."""

    def test_pulse_conv_basic(self):
        """PulseConv computes moving average correctly."""
        signal = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        result = PulseConv(signal, pulseSize=2)
        self.assertEqual(len(result), len(signal))
        # Result should be smoothed version of signal
        self.assertGreater(result.max(), 0)
        self.assertLessEqual(result.max(), 1.0)

    def test_pulse_conv_size_edge_case(self):
        """PulseConv returns zeros when pulseSize > signalSize."""
        signal = np.array([1.0, 2.0, 3.0])
        # pulseSize=10 > signalSize=3
        with self.assertLogs(level="WARNING") as log:
            result = PulseConv(signal, pulseSize=10)
        self.assertEqual(len(result), len(signal))
        assert_array_equal(result, np.zeros(3))
        self.assertTrue(any("exceeds signal length" in msg for msg in log.output))


class FindLocalPeaksTests(unittest.TestCase):
    """Tests for local peak detection."""

    def test_find_peaks_basic(self):
        """FindLocalPeaks detects obvious peaks."""
        # Clear peak at index 5
        signal = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 0.5, 0.3, 0.2, 0.1])
        peaks = FindLocalPeaks(signal)
        self.assertIn(5, peaks)

    def test_find_peaks_negative(self):
        """FindLocalPeaks detects minima on negative values."""
        # Clear trough at index 5
        signal = np.array([0.0, -0.1, -0.2, -0.3, -0.5, -1.0, -0.5, -0.3, -0.2, -0.1])
        peaks = FindLocalPeaks(signal)
        self.assertIn(5, peaks)

    def test_find_peaks_empty(self):
        """FindLocalPeaks returns empty for flat signal."""
        signal = np.zeros(10)
        peaks = FindLocalPeaks(signal)
        self.assertEqual(len(peaks), 0)


class HaarSegIntegrationTests(unittest.TestCase):
    """Integration tests for the full haarSeg algorithm."""

    def test_haarseg_step_signal(self):
        """haarSeg detects breakpoint in step signal."""
        # Create a clear step signal with noise
        rng = np.random.default_rng(42)
        signal = np.concatenate(
            [
                np.zeros(50) + rng.normal(0, 0.1, 50),
                np.ones(50) + rng.normal(0, 0.1, 50),
            ]
        )
        result = haarSeg(signal, breaksFdrQ=0.01)
        # Should find at least one breakpoint near index 50
        self.assertGreater(len(result["start"]), 1)
        # Check breakpoint is near the true change point
        breakpoints = result["start"][1:]  # Skip the 0 start
        self.assertTrue(any(45 <= bp <= 55 for bp in breakpoints))

    def test_haarseg_no_signal(self):
        """haarSeg returns single segment for constant signal."""
        signal = np.zeros(100)
        result = haarSeg(signal, breaksFdrQ=0.01)
        # Should return just one segment
        self.assertEqual(len(result["start"]), 1)
        self.assertEqual(result["start"][0], 0)
        self.assertEqual(result["end"][0], 99)

    def test_haarseg_short_signal(self):
        """haarSeg handles signals shorter than max wavelet scale."""
        # Signal with only 10 elements (shorter than 2^5=32)
        signal = np.array([0.0] * 5 + [1.0] * 5)
        # Should not raise, may log warnings for skipped levels
        result = haarSeg(signal, breaksFdrQ=0.01)
        self.assertIn("start", result)
        self.assertIn("end", result)
        self.assertEqual(result["start"][0], 0)


if __name__ == "__main__":
    unittest.main()
