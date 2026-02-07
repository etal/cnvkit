"""Supporting functions for the text/tabular-reporting commands.

Namely: breaks, genemetrics.
"""

from __future__ import annotations
import collections
import math
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import stats

from . import params, descriptives, segmetrics
from .segmetrics import segment_mean

if TYPE_CHECKING:
    from collections.abc import Iterator
    from cnvlib.cnary import CopyNumArray
    from numpy import float64


# _____________________________________________________________________________
# breaks


def do_breaks(
    probes: CopyNumArray, segments: CopyNumArray, min_probes: int = 1
) -> pd.DataFrame:
    """List the targeted genes in which a copy number breakpoint occurs.

    Parameters
    ----------
    probes : CopyNumArray
        Bin-level copy number data.
    segments : CopyNumArray
        Segmented copy number data.
    min_probes : int, optional
        Minimum number of probes required on each side of the breakpoint.
        Default is 1.

    Returns
    -------
    pd.DataFrame
        Table with columns: gene, chromosome, location, change,
        probes_left, probes_right.
    """
    intervals = get_gene_intervals(probes)
    bpoints = get_breakpoints(intervals, segments, min_probes)
    return pd.DataFrame.from_records(
        bpoints,
        columns=[
            "gene",
            "chromosome",
            "location",
            "change",
            "probes_left",
            "probes_right",
        ],
    )


def get_gene_intervals(
    all_probes: CopyNumArray, ignore: tuple[str, str, str] = params.IGNORE_GENE_NAMES
) -> collections.defaultdict[str, list[tuple[str, list[int], int]]]:
    """Tally genomic locations of each targeted gene.

    Return a dict of chromosomes to a list of tuples: (gene name, starts, end),
    where gene name is a string, starts is a sorted list of probe start
    positions, and end is the last probe's end position as an integer. (The
    endpoints are redundant since probes are adjacent.)
    """
    ignore += params.ANTITARGET_ALIASES  # type: ignore[assignment]
    # Tally the start & end points for each targeted gene; group by chromosome
    gene_probes: dict = collections.defaultdict(lambda: collections.defaultdict(list))
    for row in all_probes:
        gname = str(row.gene)
        if gname not in ignore:
            gene_probes[row.chromosome][gname].append(row)
    # Condense into a single interval for each gene
    intervals = collections.defaultdict(list)
    for chrom, gp in gene_probes.items():
        for gene, probes in gp.items():
            starts = sorted(row.start for row in probes)
            end = max(row.end for row in probes)
            intervals[chrom].append((gene, starts, end))
        intervals[chrom].sort(key=lambda gse: gse[1])
    return intervals


def get_breakpoints(
    intervals: collections.defaultdict[str, list[tuple[str, list[int], int]]],
    segments: CopyNumArray,
    min_probes: int,
) -> list[tuple[str, str, int, float64, int, int]]:
    """Identify segment breaks within the targeted intervals."""
    # TODO use segments.by_ranges(intervals)
    breakpoints = []
    for i, curr_row in enumerate(segments[:-1]):
        curr_chrom = curr_row.chromosome
        curr_end = curr_row.end
        next_row = segments[i + 1]
        # Skip if this segment is the last (or only) one on this chromosome
        if next_row.chromosome != curr_chrom:
            continue
        for gname, gstarts, gend in intervals[curr_chrom]:
            if gstarts[0] < curr_end < gend:
                probes_left = sum(s < curr_end for s in gstarts)
                probes_right = sum(s >= curr_end for s in gstarts)
                if probes_left >= min_probes and probes_right >= min_probes:
                    breakpoints.append(
                        (
                            gname,
                            curr_chrom,
                            math.ceil(curr_end),
                            next_row.log2 - curr_row.log2,
                            probes_left,
                            probes_right,
                        )
                    )
    breakpoints.sort(key=lambda row: (min(row[4], row[5]), abs(row[3])), reverse=True)
    return breakpoints


# _____________________________________________________________________________
# genemetrics


def do_genemetrics(
    cnarr: CopyNumArray,
    segments: CopyNumArray | None = None,
    threshold: float = 0.2,
    min_probes: int = 3,
    skip_low: bool = False,
    is_haploid_x_reference: bool = False,
    is_sample_female: None = None,
    diploid_parx_genome: str | None = None,
    location_stats: tuple[()] | list[str] = (),
    spread_stats: tuple[()] | list[str] = (),
    interval_stats: tuple[()] | list[str] = (),
    alpha: float = 0.05,
    bootstraps: int = 100,
    smoothed: bool | int = 10,
) -> pd.DataFrame:
    """Identify targeted genes with copy number gain or loss.

    Parameters
    ----------
    cnarr : CopyNumArray
        Bin-level copy number data.
    segments : CopyNumArray, optional
        Segmented copy number data. If provided, metrics are calculated
        per segment.
    threshold : float, optional
        Minimum absolute log2 ratio to consider a gene altered. Default is 0.2.
    min_probes : int, optional
        Minimum number of probes required to report a gene. Default is 3.
    skip_low : bool, optional
        Skip bins with low coverage. Default is False.
    is_haploid_x_reference : bool, optional
        Whether reference is male (haploid X). Default is False.
    is_sample_female : bool, optional
        Whether sample is female. If None, inferred from data.
    diploid_parx_genome : str, optional
        Reference genome name for pseudo-autosomal region handling
        (e.g., 'hg19', 'hg38', 'mm10').
    location_stats : list of str, optional
        Location statistics to compute: 'mean', 'median', 'mode', 'p_ttest'.
        Default is empty tuple.
    spread_stats : list of str, optional
        Spread statistics to compute: 'stdev', 'mad', 'mse', 'iqr', 'bivar', 'sem'.
        Default is empty tuple.
    interval_stats : list of str, optional
        Interval statistics to compute: 'ci' (confidence interval),
        'pi' (prediction interval). Default is empty tuple.
    alpha : float, optional
        Significance level for confidence/prediction intervals. Default is 0.05.
    bootstraps : int, optional
        Number of bootstrap iterations for confidence intervals. Default is 100.
    smoothed : bool or int, optional
        Smoothed bootstrap threshold for confidence intervals. If bool: True to
        always use smoothed bootstrap, False to never use it. If int: use smoothed
        bootstrap when gene has <= this many bins. Smoothed bootstrap adds
        Gaussian noise to improve CI accuracy for small genes. BCa correction
        is applied when smoothing is not used. Default is 10.

    Returns
    -------
    pd.DataFrame
        Table of genes with copy number alterations, including gene name,
        chromosome, log2 ratio, and probe counts.
    """
    if is_sample_female is None:
        is_sample_female = cnarr.guess_xx(
            is_haploid_x_reference=is_haploid_x_reference,
            diploid_parx_genome=diploid_parx_genome,
        )
    cnarr = cnarr.shift_xx(
        is_haploid_x_reference, is_sample_female, diploid_parx_genome
    )
    if segments:
        segments = segments.shift_xx(
            is_haploid_x_reference, is_sample_female, diploid_parx_genome
        )
        rows = gene_metrics_by_segment(
            cnarr,
            segments,
            threshold,
            skip_low,
            location_stats,
            spread_stats,
            interval_stats,
            alpha,
            bootstraps,
            smoothed,
        )
    else:
        rows = gene_metrics_by_gene(
            cnarr,
            threshold,
            skip_low,
            location_stats,
            spread_stats,
            interval_stats,
            alpha,
            bootstraps,
            smoothed,
        )
    rows_list = list(rows)
    columns = rows_list[0].index if len(rows_list) else cnarr._required_columns
    columns = ["gene"] + [col for col in columns if col != "gene"]
    table = pd.DataFrame.from_records(rows_list).reindex(columns=columns)
    if min_probes and len(table):
        n_probes = (
            table.segment_probes if "segment_probes" in table.columns else table.probes
        )
        table = table[n_probes >= min_probes]
    return table


def gene_metrics_by_gene(
    cnarr: CopyNumArray,
    threshold: float,
    skip_low: bool = False,
    location_stats: tuple[()] | list[str] = (),
    spread_stats: tuple[()] | list[str] = (),
    interval_stats: tuple[()] | list[str] = (),
    alpha: float = 0.05,
    bootstraps: int = 100,
    smoothed: bool | int = 10,
) -> Iterator[pd.Series]:
    """Identify genes where average bin copy ratio value exceeds `threshold`.

    NB: Adjust the sample's sex-chromosome log2 values beforehand with shift_xx,
    otherwise all chrX/chrY genes may be reported gained/lost.
    """
    for row in group_by_genes(
        cnarr,
        skip_low,
        location_stats,
        spread_stats,
        interval_stats,
        alpha,
        bootstraps,
        smoothed,
    ):
        if abs(row.log2) >= threshold and row.gene:
            yield row


def gene_metrics_by_segment(
    cnarr: CopyNumArray,
    segments: CopyNumArray,
    threshold: float,
    skip_low: bool = False,
    location_stats: tuple[()] | list[str] = (),
    spread_stats: tuple[()] | list[str] = (),
    interval_stats: tuple[()] | list[str] = (),
    alpha: float = 0.05,
    bootstraps: int = 100,
    smoothed: bool | int = 10,
) -> Iterator[pd.Series]:
    """Identify genes where segmented copy ratio exceeds `threshold`.

    In the output table, show each segment's weight and probes as segment_weight
    and segment_probes, alongside the gene-level weight and probes.

    NB: Adjust the sample's sex-chromosome log2 values beforehand with shift_xx,
    otherwise all chrX/chrY genes may be reported gained/lost.
    """
    extra_cols = [
        col
        for col in segments.data.columns
        if col not in cnarr.data.columns and col not in ("depth", "probes", "weight")
    ]
    for colname in extra_cols:
        cnarr[colname] = np.nan
    for segment, subprobes in cnarr.by_ranges(segments):
        if abs(segment.log2) >= threshold:
            for row in group_by_genes(
                subprobes,
                skip_low,
                location_stats,
                spread_stats,
                interval_stats,
                alpha,
                bootstraps,
                smoothed,
            ):
                row["log2"] = segment.log2
                if hasattr(segment, "weight"):
                    row["segment_weight"] = segment.weight
                if hasattr(segment, "probes"):
                    row["segment_probes"] = segment.probes
                for colname in extra_cols:
                    row[colname] = getattr(segment, colname)
                yield row


# ENH consolidate with CNA.squash_genes
def compute_gene_stats(
    bins: CopyNumArray,
    gene_log2: float,
    location_stats: tuple[()] | list[str] = (),
    spread_stats: tuple[()] | list[str] = (),
    interval_stats: tuple[()] | list[str] = (),
    alpha: float = 0.05,
    bootstraps: int = 100,
    smoothed: bool | int = 10,
) -> dict:
    """Compute statistics for bins within a gene.

    Similar to segmetrics.do_segmetrics, but for gene-level bins.

    Parameters
    ----------
    bins : CopyNumArray
        Bins within the gene.
    gene_log2 : float
        Gene's mean log2 value.
    location_stats : list of str, optional
        Location statistics to compute.
    spread_stats : list of str, optional
        Spread statistics to compute.
    interval_stats : list of str, optional
        Interval statistics to compute.
    alpha : float, optional
        Significance level for intervals.
    bootstraps : int, optional
        Number of bootstrap iterations.
    smoothed : bool, optional
        Use smoothed bootstrap.

    Returns
    -------
    dict
        Dictionary of computed statistics.
    """
    import warnings

    warnings.simplefilter("ignore", RuntimeWarning)

    stats_dict: dict[str, float] = {}
    if not any((location_stats, spread_stats, interval_stats)):
        return stats_dict

    bins_log2 = bins["log2"].to_numpy()
    if len(bins_log2) == 0:
        return stats_dict

    stat_funcs = {
        "mean": np.mean,
        "median": np.median,
        "mode": descriptives.modal_location,
        "p_ttest": lambda a: stats.ttest_1samp(a, 0.0, nan_policy="omit")[1],
        "stdev": np.std,
        "mad": descriptives.median_absolute_deviation,
        "mse": descriptives.mean_squared_error,
        "iqr": descriptives.interquartile_range,
        "bivar": descriptives.biweight_midvariance,
        "sem": stats.sem,
        "ci": segmetrics.make_ci_func(alpha, bootstraps, smoothed),
        "pi": segmetrics.make_pi_func(alpha),
    }

    # Location statistics
    for statname in location_stats:
        func = stat_funcs[statname]
        stats_dict[statname] = func(bins_log2)

    # Spread statistics (deviations from gene mean)
    if spread_stats:
        deviations = bins_log2 - gene_log2
        for statname in spread_stats:
            func = stat_funcs[statname]
            stats_dict[statname] = func(deviations)

    # Interval statistics
    weights = bins["weight"].to_numpy() if "weight" in bins else np.ones(len(bins_log2))
    if "ci" in interval_stats:
        ci_lo, ci_hi = stat_funcs["ci"](bins_log2, weights)
        stats_dict["ci_lo"] = ci_lo
        stats_dict["ci_hi"] = ci_hi
    if "pi" in interval_stats:
        pi_lo, pi_hi = stat_funcs["pi"](bins_log2, weights)
        stats_dict["pi_lo"] = pi_lo
        stats_dict["pi_hi"] = pi_hi

    return stats_dict


def group_by_genes(
    cnarr: CopyNumArray,
    skip_low: bool,
    location_stats: tuple[()] | list[str] = (),
    spread_stats: tuple[()] | list[str] = (),
    interval_stats: tuple[()] | list[str] = (),
    alpha: float = 0.05,
    bootstraps: int = 100,
    smoothed: bool | int = 10,
) -> Iterator[pd.Series]:
    """Group probe and coverage data by gene.

    Return an iterable of genes, in chromosomal order, associated with their
    location and coverages:

        [(gene, chrom, start, end, [coverages]), ...]
    """
    ignore = ("", np.nan, *params.ANTITARGET_ALIASES)
    for gene, rows in cnarr.by_gene():
        if not rows or gene in ignore:
            continue
        segmean = segment_mean(rows, skip_low)
        if segmean is None:
            continue
        outrow = rows[0].copy()
        outrow["end"] = rows.end.iat[-1]
        outrow["gene"] = gene
        outrow["log2"] = segmean
        outrow["probes"] = len(rows)
        if "weight" in rows:
            outrow["weight"] = rows["weight"].sum()
            if "depth" in rows:
                outrow["depth"] = np.average(rows["depth"], weights=rows["weight"])
        elif "depth" in rows:
            outrow["depth"] = rows["depth"].mean()

        # Compute statistics if requested
        if any((location_stats, spread_stats, interval_stats)):
            gene_stats = compute_gene_stats(
                rows,
                segmean,
                location_stats,
                spread_stats,
                interval_stats,
                alpha,
                bootstraps,
                smoothed,
            )
            for stat_name, stat_value in gene_stats.items():
                outrow[stat_name] = stat_value

        yield outrow
