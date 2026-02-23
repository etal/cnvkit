"""Filter copy number segments."""

from __future__ import annotations
import functools
import logging

import numpy as np
import pandas as pd

from .descriptives import weighted_median
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable
    from cnvlib.cnary import CopyNumArray
    from pandas.core.frame import DataFrame
    from pandas.core.series import Series


def require_column(*colnames) -> Callable:
    """Wrapper to coordinate the segment-filtering functions.

    Verify that the given columns are in the CopyNumArray the wrapped function
    takes. Also log the number of rows in the array before and after filtration.
    """
    if len(colnames) == 1:
        msg = "'{}' filter requires column '{}'"
    else:
        msg = "'{}' filter requires columns " + ", ".join(["'{}'"] * len(colnames))

    def wrap(func):
        @functools.wraps(func)
        def wrapped_f(segarr):
            filtname = func.__name__
            if any(c not in segarr for c in colnames):
                raise ValueError(msg.format(filtname, *colnames))
            result = func(segarr)
            logging.info(
                "Filtered by '%s' from %d to %d rows",
                filtname,
                len(segarr),
                len(result),
            )
            return result

        return wrapped_f

    return wrap


def squash_by_groups(
    cnarr: CopyNumArray, levels: Series, by_arm: bool = False
) -> CopyNumArray:
    """Reduce CopyNumArray rows to a single row within each given level."""
    # Enumerate runs of identical values
    change_levels = enumerate_changes(levels)
    assert (change_levels.index == levels.index).all()
    assert cnarr.data.index.is_unique
    assert levels.index.is_unique
    assert change_levels.index.is_unique
    if by_arm:
        # Enumerate chromosome arms
        arm_levels = []
        for i, (_chrom, cnarm) in enumerate(cnarr.by_arm()):
            arm_levels.append(np.repeat(i, len(cnarm)))
        change_levels += np.concatenate(arm_levels)
    else:
        # Enumerate chromosomes
        chrom_names = cnarr["chromosome"].unique()
        chrom_col = cnarr["chromosome"].map(
            pd.Series(np.arange(len(chrom_names)), index=chrom_names)
        )
        change_levels += chrom_col
    data = cnarr.data.assign(_group=change_levels)
    groupkey = ["_group"]
    if "cn1" in cnarr:
        # Keep allele-specific CNAs separate
        data["_g1"] = enumerate_changes(cnarr["cn1"])
        data["_g2"] = enumerate_changes(cnarr["cn2"])
        groupkey.extend(["_g1", "_g2"])
    data = (
        data.groupby(groupkey, as_index=False, group_keys=False, sort=False)[
            data.columns
        ]
        .apply(squash_region)
        .reset_index(drop=True)
    )
    return cnarr.as_dataframe(data)


def enumerate_changes(levels: Series) -> Series:
    """Assign a unique integer to each run of identical values.

    Repeated but non-consecutive values will be assigned different integers.
    """
    return levels.diff().fillna(0).abs().cumsum().astype(int)


def squash_region(cnarr: DataFrame) -> DataFrame:
    """Reduce a CopyNumArray to 1 row, keeping fields sensible.

    Most fields added by the `segmetrics` command will be dropped.
    """
    assert "weight" in cnarr
    out = {
        "chromosome": [cnarr["chromosome"].iat[0]],
        "start": cnarr["start"].iat[0],
        "end": cnarr["end"].iat[-1],
    }
    region_weight = cnarr["weight"].sum()
    if region_weight > 0:
        out["log2"] = np.average(cnarr["log2"], weights=cnarr["weight"])
    else:
        out["log2"] = np.mean(cnarr["log2"])
    out["gene"] = ",".join(cnarr["gene"].drop_duplicates())
    out["probes"] = cnarr["probes"].sum() if "probes" in cnarr else len(cnarr)
    out["weight"] = region_weight
    if "depth" in cnarr:
        if region_weight > 0:
            out["depth"] = np.average(cnarr["depth"], weights=cnarr["weight"])
        else:
            out["depth"] = np.mean(cnarr["depth"])
    if "baf" in cnarr:
        if region_weight > 0:
            out["baf"] = np.average(cnarr["baf"], weights=cnarr["weight"])
        else:
            out["baf"] = np.mean(cnarr["baf"])
    if "cn" in cnarr:
        if region_weight > 0:
            out["cn"] = weighted_median(cnarr["cn"], cnarr["weight"])
        else:
            out["cn"] = np.median(cnarr["cn"])
        if "cn1" in cnarr:
            if region_weight > 0:
                out["cn1"] = weighted_median(cnarr["cn1"], cnarr["weight"])
            else:
                out["cn1"] = np.median(cnarr["cn1"])
            out["cn2"] = out["cn"] - out["cn1"]
    if "p_bintest" in cnarr:
        # Only relevant for single-bin segments, but this seems safe/conservative
        out["p_bintest"] = cnarr["p_bintest"].max()
    return pd.DataFrame(out)


@require_column("cn")
def ampdel(segarr: CopyNumArray) -> CopyNumArray:
    """Merge segments by amplified/deleted/neutral copy number status.

    Follow the clinical reporting convention:

    - 5+ copies (2.5-fold gain) is amplification
    - 0 copies is homozygous/deep deletion
    - CNAs of lesser degree are not reported

    This is recommended only for selecting segments corresponding to
    actionable, usually focal, CNAs. Any real and potentially informative but
    lower-level CNAs will be dropped.
    """
    levels = np.zeros(len(segarr))
    levels[segarr["cn"] == 0] = -1
    levels[segarr["cn"] >= 5] = 1
    # or: segarr['log2'] >= np.log2(2.5)
    cnarr = squash_by_groups(segarr, pd.Series(levels))
    return cnarr[(cnarr["cn"] == 0) | (cnarr["cn"] >= 5)]


@require_column("probes")
def bic(segarr: CopyNumArray) -> CopyNumArray:
    """Merge adjacent segments whose difference is not justified by BIC.

    Uses the Bayesian Information Criterion to test whether two adjacent
    same-chromosome segments are better modeled as separate or merged.
    Iteratively merges the pair with the most negative delta-BIC until
    no pair benefits from merging.

    Variance is estimated from the ``stdev`` column if present (from
    ``segmetrics --spread stdev``), otherwise from the ``weight`` column.

    Complexity is O(k^2) in the number of segments per chromosome in the worst
    case (each iteration scans all adjacent pairs). This is fast for typical
    segment counts (tens to low hundreds) but may be slow for very heavily
    segmented inputs.

    See: BIC-seq (Xi 2011), doi:10.1073/pnas.1110574108
    """
    if len(segarr) < 2:
        return segarr

    data = segarr.data

    # --- Determine per-segment RSS values ---
    n = data["probes"].to_numpy(dtype=float)
    mu = data["log2"].to_numpy(dtype=float)

    if "stdev" in data.columns:
        sd = data["stdev"].to_numpy(dtype=float).copy()
        # Fallback for segments with stdev=0 or probes=1: use median stdev
        needs_fallback = (sd == 0) | (n <= 1)
        valid = ~needs_fallback
        if valid.any():
            fallback_sd = float(np.median(sd[valid]))
        else:
            fallback_sd = 0.0
        sd[needs_fallback] = fallback_sd
        # stdev is the SD of individual bin log2 values within the segment,
        # so RSS = variance * n_observations = stdev^2 * probes
        rss = sd**2 * n
    elif "weight" in data.columns:
        w = data["weight"].to_numpy(dtype=float)
        # Variance ~ n/weight, so RSS = n * variance = n^2/weight
        needs_fallback = w == 0
        valid = ~needs_fallback
        if valid.any():
            fallback_rss = float(np.median(n[valid] ** 2 / w[valid]))
        else:
            fallback_rss = 0.0
        rss = np.where(needs_fallback, fallback_rss, n**2 / w)
    else:
        logging.warning("BIC filter: no 'stdev' or 'weight' column; skipping")
        return segarr

    chrom = data["chromosome"].to_numpy()

    # --- Iterative greedy merging ---
    # Track groups of original row indices, plus stats for BIC computation
    groups: list[list[int]] = [[j] for j in range(len(data))]
    n_list = list(n)
    mu_list = list(mu)
    rss_list = list(rss)
    chrom_list = list(chrom)

    def _delta_bic(i: int) -> float:
        """Compute delta-BIC for merging segments at positions i and i+1."""
        n1, n2 = n_list[i], n_list[i + 1]
        mu1, mu2 = mu_list[i], mu_list[i + 1]
        rss1, rss2 = rss_list[i], rss_list[i + 1]

        n_tot = n1 + n2
        rss_sep = rss1 + rss2

        # Additional RSS from merging means
        rss_add = n1 * n2 / n_tot * (mu1 - mu2) ** 2
        rss_merged = rss_sep + rss_add

        if rss_sep == 0:
            # Both segments have zero within-segment variance
            if mu1 == mu2:
                return float(-np.log(n_tot))  # Merge: saves a parameter
            return np.inf  # Different means with zero noise: don't merge
        return float(n_tot * np.log(rss_merged / rss_sep) - np.log(n_tot))

    while True:
        # Find the pair with the most negative delta-BIC (same chromosome)
        best_dbic = 0.0
        best_idx = -1
        for i in range(len(n_list) - 1):
            if chrom_list[i] != chrom_list[i + 1]:
                continue
            dbic = _delta_bic(i)
            if dbic < best_dbic:
                best_dbic = dbic
                best_idx = i

        if best_idx < 0:
            break  # No beneficial merges remain

        # Merge groups[best_idx] and groups[best_idx + 1]
        i = best_idx
        new_n = n_list[i] + n_list[i + 1]
        new_mu = (n_list[i] * mu_list[i] + n_list[i + 1] * mu_list[i + 1]) / new_n
        rss_add = n_list[i] * n_list[i + 1] / new_n * (mu_list[i] - mu_list[i + 1]) ** 2
        new_rss = rss_list[i] + rss_list[i + 1] + rss_add

        groups[i] = groups[i] + groups[i + 1]
        n_list[i] = new_n
        mu_list[i] = new_mu
        rss_list[i] = new_rss

        del groups[i + 1]
        del n_list[i + 1]
        del mu_list[i + 1]
        del rss_list[i + 1]
        del chrom_list[i + 1]

    # Build result by squashing each group of original rows
    result_rows = []
    for group in groups:
        result_rows.append(squash_region(data.iloc[group]))
    result = pd.concat(result_rows, ignore_index=True)
    return segarr.as_dataframe(result)


@require_column("ci_lo", "ci_hi")
def ci(segarr: CopyNumArray) -> CopyNumArray:
    """Merge segments by confidence interval (overlapping 0).

    Segments with lower CI above 0 are kept as gains, upper CI below 0 as
    losses, and the rest with CI overlapping zero are collapsed as neutral.

    Only adjacent neutral segments (CI overlapping zero) are merged together.
    Segments that are confidently non-neutral are each preserved individually,
    even if adjacent segments have the same sign.
    """
    levels = np.zeros(len(segarr))
    levels[segarr["ci_lo"].to_numpy() > 0] = 1
    levels[segarr["ci_hi"].to_numpy() < 0] = -1
    # Assign each non-neutral segment a unique level so that squash_by_groups
    # never merges two non-neutral segments together, even if adjacent.  Without
    # this, two adjacent losses (both level -1) would be squashed into one,
    # losing the distinct log2 magnitudes (e.g. -1.0 and -0.03).
    nonzero_mask = levels != 0
    if nonzero_mask.any():
        levels[nonzero_mask] = np.arange(1, nonzero_mask.sum() + 1)
    return squash_by_groups(segarr, pd.Series(levels))


@require_column("cn")
def cn(segarr: CopyNumArray) -> CopyNumArray:
    """Merge segments by integer copy number."""
    return squash_by_groups(segarr, segarr["cn"])


@require_column("sem")
def sem(segarr: CopyNumArray, zscore: float = 1.96) -> CopyNumArray:
    """Merge segments by Standard Error of the Mean (SEM).

    Use each segment's SEM value to estimate a 95% confidence interval (via
    `zscore`). Segments with lower CI above 0 are kept as gains, upper CI below
    0 as losses, and the rest with CI overlapping zero are collapsed as neutral.

    Only adjacent neutral segments (CI overlapping zero) are merged together.
    Segments that are confidently non-neutral are each preserved individually,
    even if adjacent segments have the same sign.
    """
    margin = segarr["sem"] * zscore
    levels = np.zeros(len(segarr))
    levels[segarr["log2"] - margin > 0] = 1
    levels[segarr["log2"] + margin < 0] = -1
    # Assign each non-neutral segment a unique level so that squash_by_groups
    # never merges two non-neutral segments together, even if adjacent
    nonzero_mask = levels != 0
    if nonzero_mask.any():
        levels[nonzero_mask] = np.arange(1, nonzero_mask.sum() + 1)
    return squash_by_groups(segarr, pd.Series(levels))
