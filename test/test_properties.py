"""Property-based tests for CNVkit numerical invariants using Hypothesis.

Tests mathematical properties (equivariance, boundedness, round-trips) for:
- cnvlib/descriptives.py  -- robust statistical estimators
- cnvlib/smoothing.py     -- signal-smoothing functions
- cnvlib/call.py          -- copy number arithmetic helpers

See: https://github.com/etal/cnvkit/issues/1038
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from cnvlib.call import (
    _log2_ratio_to_absolute,
    _log2_ratio_to_absolute_pure,
    rescale_baf,
)
from cnvlib.descriptives import (
    biweight_location,
    biweight_midvariance,
    gapper_scale,
    interquartile_range,
    median_absolute_deviation,
    q_n,
    weighted_median,
)
from cnvlib.smoothing import (
    _width2wing,
    kaiser,
    rolling_median,
    savgol,
)

# ---------------------------------------------------------------------------
# Hypothesis settings profiles
# ---------------------------------------------------------------------------

settings.register_profile(
    "ci",
    max_examples=50,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.differing_executors],
    deadline=500,
)
settings.register_profile(
    "default",
    max_examples=100,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.differing_executors],
    deadline=800,
)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))

pytestmark = pytest.mark.hypothesis


# ---------------------------------------------------------------------------
# Shared strategy helpers
# ---------------------------------------------------------------------------

# Finite floats in the plausible CNV log2 ratio range
_LOG2_FLOATS = st.floats(
    min_value=-8.0, max_value=8.0, allow_nan=False, allow_infinity=False
)

# Positive finite floats (for scale factors, ref_copies)
_POSITIVE_FLOATS = st.floats(
    min_value=1e-3, max_value=1e4, allow_nan=False, allow_infinity=False
)


def _float_array(min_size: int = 2, max_size: int = 200):
    """Strategy: list of finite log2-range floats, converted to np.array."""
    return st.lists(_LOG2_FLOATS, min_size=min_size, max_size=max_size).map(np.array)


def _non_constant_array(min_size: int = 3, max_size: int = 200):
    """Strategy: finite float array with at least two distinct values."""
    return _float_array(min_size, max_size).filter(lambda a: np.ptp(a) > 1e-12)


def _weighted_pair(min_size: int = 2, max_size: int = 100):
    """Strategy: paired (values, weights) arrays of equal length."""
    return st.integers(min_value=min_size, max_value=max_size).flatmap(
        lambda n: st.tuples(
            st.lists(_LOG2_FLOATS, min_size=n, max_size=n).map(np.array),
            st.lists(
                st.floats(
                    min_value=0.01,
                    max_value=1e3,
                    allow_nan=False,
                    allow_infinity=False,
                ),
                min_size=n,
                max_size=n,
            ).map(np.array),
        )
    )


def _smoothing_input(min_n: int = 10, max_n: int = 100):
    """Strategy: (array, fractional_width) for smoothing functions."""
    return st.integers(min_value=min_n, max_value=max_n).flatmap(
        lambda n: st.tuples(
            st.lists(_LOG2_FLOATS, min_size=n, max_size=n).map(np.array),
            st.floats(
                min_value=0.05, max_value=0.9, allow_nan=False, allow_infinity=False
            ),
        )
    )


def _smoothing_input_int(min_n: int = 10, max_n: int = 100):
    """Strategy: (array, integer_width) for smoothing functions."""
    return st.integers(min_value=min_n, max_value=max_n).flatmap(
        lambda n: st.tuples(
            st.lists(_LOG2_FLOATS, min_size=n, max_size=n).map(np.array),
            st.integers(min_value=4, max_value=max(4, n - 1)),
        )
    )


# ---------------------------------------------------------------------------
# Tests: cnvlib/descriptives.py
# ---------------------------------------------------------------------------


class TestBiweightLocation:
    """Invariants for biweight_location (robust location estimator)."""

    @given(_float_array())
    def test_finite_output(self, a):
        result = biweight_location(a)
        assert np.isfinite(result)

    @given(_float_array())
    def test_bounded_by_data_range(self, a):
        result = biweight_location(a)
        # Allow tiny float rounding (iterative estimator can exceed range by ~eps)
        eps = 1e-12
        assert a.min() - eps <= result <= a.max() + eps

    @given(_non_constant_array(), _LOG2_FLOATS)
    def test_shift_equivariance(self, a, c):
        """biweight_location(a + c) == biweight_location(a) + c."""
        assume(np.all(np.isfinite(a + c)))
        base = biweight_location(a)
        shifted = biweight_location(a + c)
        tol = 1e-4 * (abs(base) + abs(c) + 1)
        assert abs(shifted - (base + c)) < tol

    @given(
        # min_size=20: biweight outlier-rejection weights are unstable for
        # small samples, breaking scale equivariance even for well-conditioned
        # data.  20+ points give enough mass for stable convergence.
        _non_constant_array(min_size=20),
        st.floats(min_value=0.5, max_value=10.0, allow_nan=False, allow_infinity=False),
    )
    def test_scale_equivariance(self, a, k):
        """biweight_location(a * k) == biweight_location(a) * k."""
        assume(np.all(np.isfinite(a * k)))
        base = biweight_location(a)
        scaled = biweight_location(a * k)
        # Iterative convergence (epsilon=1e-3 per run) + floating-point
        # arithmetic across scales can compound to ~3 * epsilon.
        tol = max(0.02 * abs(base * k), 0.004)
        assert abs(scaled - base * k) < tol


class TestBiweightMidvariance:
    """Invariants for biweight_midvariance (robust scale estimator)."""

    @given(_float_array())
    def test_non_negative(self, a):
        assert biweight_midvariance(a) >= 0

    @given(_LOG2_FLOATS, st.integers(min_value=2, max_value=50))
    def test_zero_for_constant_array(self, c, n):
        a = np.full(n, c)
        assert biweight_midvariance(a) == 0

    @given(
        _non_constant_array(min_size=20),
        st.floats(min_value=0.5, max_value=10.0, allow_nan=False, allow_infinity=False),
    )
    def test_scale_equivariance(self, a, k):
        """biweight_midvariance(a * k) == biweight_midvariance(a) * k."""
        assume(np.all(np.isfinite(a * k)))
        base = biweight_midvariance(a)
        scaled = biweight_midvariance(a * k)
        tol = max(0.02 * abs(base * k), 0.004)
        assert abs(scaled - base * k) < tol


class TestMedianAbsoluteDeviation:
    """Invariants for median_absolute_deviation."""

    @given(_float_array())
    def test_non_negative(self, a):
        assert median_absolute_deviation(a) >= 0

    @given(_float_array(), _LOG2_FLOATS)
    def test_shift_invariance(self, a, c):
        """MAD is location-free: mad(a + c) == mad(a)."""
        assume(np.all(np.isfinite(a + c)))
        assert (
            abs(median_absolute_deviation(a + c) - median_absolute_deviation(a)) < 1e-10
        )

    @given(_non_constant_array(), _POSITIVE_FLOATS)
    def test_scale_equivariance(self, a, k):
        assume(np.all(np.isfinite(a * k)))
        base = median_absolute_deviation(a)
        scaled = median_absolute_deviation(a * k)
        assert abs(scaled - base * k) < 1e-8 * (base * k + 1)

    @given(_float_array())
    def test_scale_to_sd_relationship(self, a):
        """scale_to_sd=True equals scale_to_sd=False * 1.4826."""
        raw = median_absolute_deviation(a, scale_to_sd=False)
        scaled = median_absolute_deviation(a, scale_to_sd=True)
        assert abs(scaled - raw * 1.4826) < 1e-10


class TestInterquartileRange:
    """Invariants for interquartile_range."""

    @given(_float_array())
    def test_non_negative(self, a):
        assert interquartile_range(a) >= 0

    @given(_float_array(), _LOG2_FLOATS)
    def test_shift_invariance(self, a, c):
        assume(np.all(np.isfinite(a + c)))
        assert abs(interquartile_range(a + c) - interquartile_range(a)) < 1e-10

    @given(_non_constant_array(), _POSITIVE_FLOATS)
    def test_scale_equivariance(self, a, k):
        assume(np.all(np.isfinite(a * k)))
        base = interquartile_range(a)
        scaled = interquartile_range(a * k)
        assert abs(scaled - base * k) < 1e-8 * (base * k + 1)

    @given(_LOG2_FLOATS, st.integers(min_value=2, max_value=50))
    def test_zero_for_constant_array(self, c, n):
        assert interquartile_range(np.full(n, c)) == 0


class TestWeightedMedian:
    """Invariants for weighted_median."""

    @given(_weighted_pair())
    def test_bounded_by_data_range(self, pair):
        a, w = pair
        assume(w.sum() > 0)
        result = weighted_median(a, w)
        assert a.min() <= result <= a.max()

    @given(_float_array(min_size=3))
    def test_uniform_weights_match_median(self, a):
        """Uniform weights produce the same result as np.median."""
        w = np.ones(len(a))
        result = weighted_median(a, w)
        expected = np.median(a)
        assert abs(result - expected) < 1e-10

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError, match="Unequal array lengths"):
            weighted_median(np.array([1.0, 2.0]), np.array([1.0]))


class TestGapperScale:
    """Invariants for gapper_scale."""

    @given(_float_array())
    def test_non_negative(self, a):
        assert gapper_scale(a) >= 0

    @given(_LOG2_FLOATS, st.integers(min_value=2, max_value=20))
    def test_zero_for_constant_array(self, c, n):
        assert gapper_scale(np.full(n, c)) == 0

    @given(_non_constant_array(), _POSITIVE_FLOATS)
    def test_scale_equivariance(self, a, k):
        assume(np.all(np.isfinite(a * k)))
        base = gapper_scale(a)
        scaled = gapper_scale(a * k)
        tol = 1e-6 * (base * k + 1)
        assert abs(scaled - base * k) < tol


class TestQn:
    """Invariants for the Rousseeuw-Croux Q_n scale estimator."""

    @given(_float_array(max_size=50))
    def test_non_negative(self, a):
        assert q_n(a) >= 0

    @given(_non_constant_array(max_size=40), _POSITIVE_FLOATS)
    def test_scale_equivariance(self, a, k):
        assume(np.all(np.isfinite(a * k)))
        base = q_n(a)
        scaled = q_n(a * k)
        tol = 1e-6 * (base * k + 1)
        assert abs(scaled - base * k) < tol


# ---------------------------------------------------------------------------
# Tests: cnvlib/smoothing.py
# ---------------------------------------------------------------------------


class TestWidth2Wing:
    """Invariants for _width2wing (fractional/integer width to half-window)."""

    @given(
        st.integers(min_value=10, max_value=200),
        st.floats(
            min_value=0.01, max_value=0.99, allow_nan=False, allow_infinity=False
        ),
    )
    def test_fractional_width_returns_bounded_int(self, n, frac):
        x = np.zeros(n)
        wing = _width2wing(frac, x)
        assert isinstance(wing, int)
        assert 1 <= wing <= n - 1

    @given(st.integers(min_value=10, max_value=200))
    def test_integer_width_returns_bounded_int(self, n):
        x = np.zeros(n)
        wing = _width2wing(4, x)
        assert isinstance(wing, int)
        assert 1 <= wing <= n - 1

    @given(st.integers(min_value=10, max_value=200))
    def test_invalid_width_raises(self, n):
        x = np.zeros(n)
        with pytest.raises(ValueError, match="fraction between 0 and 1"):
            _width2wing(1.5, x)


class TestRollingMedian:
    """Invariants for rolling_median."""

    @given(_smoothing_input())
    def test_output_length_preserved(self, args):
        x, width = args
        result = rolling_median(x, width)
        assert len(result) == len(x)

    @given(_smoothing_input())
    def test_output_bounded_by_input_range(self, args):
        x, width = args
        result = rolling_median(x, width)
        assert result.min() >= x.min() - 1e-10
        assert result.max() <= x.max() + 1e-10

    @given(_smoothing_input())
    def test_no_nans_in_output(self, args):
        x, width = args
        result = rolling_median(x, width)
        assert not np.any(np.isnan(result))


class TestKaiser:
    """Invariants for kaiser smoothing."""

    @given(_smoothing_input())
    def test_output_length_preserved(self, args):
        x, width = args
        result = kaiser(x, width)
        assert len(result) == len(x)

    @given(_smoothing_input())
    def test_no_nans_in_output(self, args):
        x, width = args
        result = kaiser(x, width)
        assert not np.any(np.isnan(result))


class TestSavgol:
    """Invariants for Savitzky-Golay smoothing."""

    @given(_smoothing_input_int())
    def test_output_length_preserved(self, args):
        x, width = args
        result = savgol(x, total_width=width)
        assert len(result) == len(x)

    @given(
        st.integers(min_value=10, max_value=100),
        _LOG2_FLOATS,
        st.integers(min_value=4, max_value=20),
    )
    def test_constant_array_unchanged(self, n, c, width):
        """Savgol of a constant array returns the same constant."""
        width = min(width, n - 1)
        x = np.full(n, c)
        result = savgol(x, total_width=width)
        assert np.allclose(result, c, atol=1e-6)


# ---------------------------------------------------------------------------
# Tests: cnvlib/call.py (pure arithmetic helpers)
# ---------------------------------------------------------------------------


class TestLog2RatioToAbsolutePure:
    """Invariants for _log2_ratio_to_absolute_pure."""

    @given(_LOG2_FLOATS, _POSITIVE_FLOATS)
    def test_positive_output(self, log2_ratio, ref_copies):
        result = _log2_ratio_to_absolute_pure(log2_ratio, ref_copies)
        assert result > 0

    @given(_POSITIVE_FLOATS)
    def test_zero_log2_returns_ref_copies(self, ref_copies):
        """log2=0 means no change: result == ref_copies."""
        result = _log2_ratio_to_absolute_pure(0.0, ref_copies)
        assert abs(result - ref_copies) < 1e-12

    @given(_LOG2_FLOATS, _POSITIVE_FLOATS)
    def test_round_trip(self, log2_ratio, ref_copies):
        """np.log2(result / ref_copies) recovers the log2 ratio."""
        result = _log2_ratio_to_absolute_pure(log2_ratio, ref_copies)
        recovered = np.log2(result / ref_copies)
        assert abs(recovered - log2_ratio) < 1e-9

    @given(
        st.tuples(_LOG2_FLOATS, _LOG2_FLOATS).filter(
            lambda pair: pair[1] - pair[0] > 1e-10
        ),
        _POSITIVE_FLOATS,
    )
    def test_monotone_in_log2(self, pair, ref_copies):
        """Larger log2 ratio always gives more copies."""
        v1, v2 = pair
        r1 = _log2_ratio_to_absolute_pure(v1, ref_copies)
        r2 = _log2_ratio_to_absolute_pure(v2, ref_copies)
        assert r1 < r2


class TestLog2RatioToAbsolute:
    """Invariants for _log2_ratio_to_absolute (purity-corrected)."""

    @given(_LOG2_FLOATS, _POSITIVE_FLOATS, _POSITIVE_FLOATS)
    def test_purity_one_equals_pure(self, log2_ratio, ref_copies, expect_copies):
        """With purity=1, the impure formula reduces to the pure formula."""
        impure = _log2_ratio_to_absolute(log2_ratio, ref_copies, expect_copies, 1.0)
        pure = _log2_ratio_to_absolute_pure(log2_ratio, ref_copies)
        assert abs(impure - pure) < 1e-9

    @given(
        _POSITIVE_FLOATS,
        st.floats(
            min_value=0.05, max_value=0.99, allow_nan=False, allow_infinity=False
        ),
    )
    def test_neutral_log2_gives_expected_copies(self, copies, purity):
        """When log2=0, autosomal result equals the expected copy number."""
        # For autosomes: ref_copies == expect_copies, neutral log2 = 0
        result = _log2_ratio_to_absolute(0.0, copies, copies, purity)
        assert abs(result - copies) < 1e-9


class TestRescaleBaf:
    """Invariants for rescale_baf (purity-adjusted B-allele frequency)."""

    @given(
        st.lists(
            st.floats(
                min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False
            ),
            min_size=1,
            max_size=50,
        ).map(pd.Series),
    )
    def test_purity_one_is_identity(self, baf):
        """rescale_baf with purity=1 returns the input unchanged."""
        result = rescale_baf(1.0, baf)
        assert np.allclose(result.values, baf.values, atol=1e-12)

    @given(
        st.floats(
            min_value=0.05, max_value=0.99, allow_nan=False, allow_infinity=False
        ),
        st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
    )
    def test_mixture_round_trip(self, purity, tumor_baf):
        """Mixing tumor+normal BAF, then rescaling, recovers tumor BAF."""
        observed = purity * tumor_baf + (1 - purity) * 0.5
        recovered = rescale_baf(purity, pd.Series([observed]))
        assert abs(recovered.iloc[0] - tumor_baf) < 1e-9

    @given(
        st.floats(
            min_value=0.05, max_value=0.99, allow_nan=False, allow_infinity=False
        ),
    )
    def test_balanced_baf_stays_at_half(self, purity):
        """Observed BAF=0.5 always rescales to 0.5 (balanced heterozygosity)."""
        result = rescale_baf(purity, pd.Series([0.5]))
        assert abs(result.iloc[0] - 0.5) < 1e-12

    @given(
        st.floats(
            min_value=0.05, max_value=0.99, allow_nan=False, allow_infinity=False
        ),
        st.lists(
            st.floats(
                min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False
            ),
            min_size=1,
            max_size=50,
        ).map(pd.Series),
    )
    def test_output_clamped_to_unit_interval(self, purity, observed_baf):
        """rescale_baf output is always in [0, 1] for any in-range inputs."""
        result = rescale_baf(purity, observed_baf)
        assert (result.values >= 0.0).all()
        assert (result.values <= 1.0).all()

    def test_low_purity_extreme_baf_clamped_to_zero(self):
        """Regression for issue #601: low purity + observed BAF ~ 0 must
        clamp to 0, not produce a negative value.
        """
        # purity=0.37 (from issue #601). Without clamping:
        # tumor_baf = (0.0 - 0.5*0.63) / 0.37 = -0.851
        result = rescale_baf(0.37, pd.Series([0.0]))
        assert result.iloc[0] == 0.0

    def test_low_purity_extreme_baf_clamped_to_one(self):
        """Symmetric case: observed BAF ~ 1 must clamp to 1, not exceed it."""
        # tumor_baf = (1.0 - 0.5*0.63) / 0.37 = 1.851
        result = rescale_baf(0.37, pd.Series([1.0]))
        assert result.iloc[0] == 1.0

    def test_nan_preserved(self):
        """NaN BAF (segments with no SNP coverage) passes through unchanged."""
        result = rescale_baf(0.5, pd.Series([0.4, np.nan, 0.6]))
        assert not np.isnan(result.iloc[0])
        assert np.isnan(result.iloc[1])
        assert not np.isnan(result.iloc[2])

    def test_no_warning_at_float_noise_boundary(self, caplog):
        """Values within float-arithmetic noise of 0 or 1 must not log a warning."""
        # observed_baf at the exact lower boundary 0.5*(1-purity) gives tumor_baf=0
        purity = 0.5
        boundary_low = 0.5 * (1 - purity)
        boundary_high = 0.5 + 0.5 * purity
        with caplog.at_level("WARNING"):
            rescale_baf(purity, pd.Series([boundary_low, 0.5, boundary_high]))
        assert not any("clamped" in rec.message for rec in caplog.records)
