"""Segmentation of copy number values."""

from __future__ import annotations

import locale
import logging
import tempfile
from io import StringIO
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from skgenome import tabio
from skgenome.combiners import join_strings
from skgenome.intersect import iter_slices

from .. import core, parallel, params, smoothing
from ..cnary import CopyNumArray as CNA
from . import cbs, haar, hmm, none

if TYPE_CHECKING:
    from cnvlib.vary import VariantArray


SEGMENT_METHODS = ("cbs", "haar", "none", "hmm", "hmm-tumor", "hmm-germline")
# HMM method variants, for dispatch (replaces scattered method.startswith("hmm"))
_HMM_METHODS = frozenset({"hmm", "hmm-tumor", "hmm-germline"})


def do_segmentation(
    cnarr: CNA,
    method: str,
    diploid_parx_genome: str | None = None,
    threshold: float | None = None,
    variants: VariantArray | None = None,
    skip_low: bool = False,
    skip_outliers: int = 10,
    min_weight: int = 0,
    save_dataframe: bool = False,
    rscript_path: str = "Rscript",
    processes: int = 1,
    smooth_cbs: bool = False,
) -> CNA | tuple[CNA, str]:
    """Infer copy number segments from the given coverage table.

    Parameters
    ----------
    cnarr : CopyNumArray
        Bin-level copy number ratios (.cnr file).
    method : str
        Segmentation algorithm: 'cbs', 'haar', 'hmm', 'hmm-tumor',
        or 'hmm-germline'.
    diploid_parx_genome : str, optional
        Reference genome name for pseudo-autosomal region handling
        (e.g., 'hg19', 'hg38', 'mm10').
    threshold : float, optional
        Significance threshold (for CBS/haar) or smoothing window size
        (for HMM methods). If None, uses method-specific defaults.
    variants : VariantArray, optional
        Variant allele frequencies to incorporate into HMM segmentation.
    skip_low : bool, optional
        Skip bins with low coverage. Default is False.
    skip_outliers : int, optional
        Skip bins with log2 ratios more than this many standard deviations
        from the chromosome arm mean. Default is 10.
    min_weight : int, optional
        Minimum weight threshold for including bins. Default is 0.
    save_dataframe : bool, optional
        Return the R dataframe as a string along with segments. Default is False.
    rscript_path : str, optional
        Path to Rscript executable for the CBS method. Default is "Rscript".
    processes : int, optional
        Number of parallel processes to use. Default is 1.
    smooth_cbs : bool, optional
        Apply smoothing to CBS segmentation results. Default is False.

    Returns
    -------
    CopyNumArray or tuple
        Segmented copy number data (.cns format). If save_dataframe=True,
        returns (segments, R dataframe string).

    Raises
    ------
    ValueError
        If method is not one of the supported segmentation methods.
    """
    if method not in SEGMENT_METHODS:
        raise ValueError(
            "'method' must be one of: "
            + ", ".join(SEGMENT_METHODS)
            + "; got: "
            + repr(method)
        )

    # Handle empty input (e.g., CNR file with only header row)
    if not len(cnarr):
        if save_dataframe:
            return cnarr, ""
        return cnarr

    if cnarr.sample_id is None:
        logging.warning(
            "Input has no sample_id set; the segmented output will be "
            "unlabeled (CLI usage derives sample_id from the input filename)."
        )

    if variants is not None and len(variants):
        # Process variants in genomic order. by_ranges/baf_by_ranges slice via
        # binary search and the HMM expects an ordered observation sequence, so
        # an unsorted input VCF -- which nothing upstream sorts -- mis-assigns
        # variants to segments (silently wrong BAF) and crashes downstream in
        # variants_in_segment with "Improper post-processing of segment"
        # (#893). Copy first so we don't mutate the caller's array.
        variants = variants.copy()
        variants.sort()

    if threshold is None:
        threshold = {
            "cbs": 0.0001,
            "haar": 0.0001,
        }.get(method)
    msg = "Segmenting with method " + repr(method)
    if threshold is not None:
        if method in _HMM_METHODS:
            msg += f", smoothing window size {threshold},"
        else:
            msg += f", significance threshold {threshold},"
    msg += f" in {processes} processes"
    logging.info(msg)

    # HMM methods fit a single model across all bins (whole-genome Viterbi
    # decode and parameter estimation), so they run on the full array rather
    # than per chromosome arm in parallel.
    if method in _HMM_METHODS:
        # ENH segment p/q arms separately
        # -> assign separate identifiers via chrom name suffix?
        cna = _do_segmentation(
            cnarr,
            method,
            diploid_parx_genome,
            threshold,
            variants,
            skip_low,
            skip_outliers,
            min_weight,
            save_dataframe,
            rscript_path,
            smooth_cbs,
        )
        if save_dataframe:
            cna, rstr = cna  # type: ignore[misc]
            rstr = _to_str(rstr)

    else:
        with parallel.pick_pool(processes) as pool:
            rets = list(
                pool.map(
                    _ds,
                    (
                        (
                            ca,
                            method,
                            diploid_parx_genome,
                            threshold,
                            variants,
                            skip_low,
                            skip_outliers,
                            min_weight,
                            save_dataframe,
                            rscript_path,
                            smooth_cbs,
                        )
                        for _, ca in cnarr.by_arm()
                    ),
                )
            )
        if save_dataframe:
            # rets is a list of (CNA, R dataframe string) -- unpack
            rets, r_dframe_strings = zip(*rets, strict=True)  # type: ignore[assignment]
            # Strip the header line from all but the first dataframe, then combine
            r_dframe_iter = map(_to_str, r_dframe_strings)
            rstr = [next(r_dframe_iter)]  # type: ignore[arg-type]
            rstr.extend(r[r.index("\n") + 1 :] for r in r_dframe_iter)
            rstr = "".join(rstr)
        cna = cnarr.concat(rets)

    cna.sort_columns()  # type: ignore[union-attr]
    if save_dataframe:
        return cna, rstr  # type: ignore[return-value]
    return cna


def _to_str(s, enc=locale.getpreferredencoding()):  # noqa: B008
    # Evaluate encoding once at function definition for performance
    if isinstance(s, bytes):
        return s.decode(enc)
    return s


def _ds(
    args: tuple[
        CNA,
        str,
        str | None,
        float | None,
        VariantArray | None,
        bool,
        int,
        int,
        bool,
        str,
        bool,
    ],
) -> CNA | tuple[CNA, str]:
    """Wrapper for parallel map"""
    return _do_segmentation(*args)


def _do_segmentation(
    cnarr: CNA,
    method: str,
    diploid_parx_genome: str | None,
    threshold: float | None,
    variants: VariantArray | None = None,
    skip_low: bool = False,
    skip_outliers: int = 10,
    min_weight: int = 0,
    save_dataframe: bool = False,
    rscript_path: str = "Rscript",
    smooth_cbs: bool = False,
) -> CNA | tuple[CNA, str]:
    """Infer copy number segments from the given coverage table."""
    if not len(cnarr):
        return cnarr

    filtered_cn = cnarr.copy()
    # Drop bins with non-finite log2 (NaN or +/-inf) before any analysis
    # (#881, #508). They carry no signal and cannot be segmented, but a
    # bare comparison treats NaN as False, so without this they survive the
    # gates below and reach either the Savitzky-Golay outlier filter (scipy's
    # lstsq rejects non-finite input with "array must not contain infs or
    # NaNs") or DNAcopy's CBS -- both of which then crash. ``.isna()`` only
    # catches NaN; broadening to ``~np.isfinite`` covers the flat-reference
    # WGS path of #508 where degenerate reference subtraction can yield
    # +/-inf. read_tab drops them on the file path; do it here too for the
    # in-memory/API path (e.g. batch).
    log2_bad = ~np.isfinite(filtered_cn["log2"])
    n_log2_bad = int(log2_bad.sum())
    if n_log2_bad:
        filtered_cn = filtered_cn[~log2_bad]
        logging.info("Dropped %d bins with non-finite log2 values", n_log2_bad)
    # Filter out bins with no or near-zero sequencing coverage
    if skip_low:
        filtered_cn = filtered_cn.drop_low_coverage(verbose=False)
    # Filter by distance from rolling quantiles
    if skip_outliers:
        filtered_cn = drop_outliers(filtered_cn, 50, skip_outliers)
    # Filter by bin weights
    if min_weight:
        weight_too_low = (filtered_cn["weight"] < min_weight).fillna(True)
    else:
        weight_too_low = (filtered_cn["weight"] == 0).fillna(True)
    n_weight_too_low = weight_too_low.sum() if len(weight_too_low) else 0
    if n_weight_too_low:
        filtered_cn = filtered_cn[~weight_too_low]
        if min_weight:
            logging.debug(
                "Dropped %d bins with weight below %s", n_weight_too_low, min_weight
            )
        else:
            logging.debug("Dropped %d bins with zero weight", n_weight_too_low)

    if len(filtered_cn) != len(cnarr):
        msg = f"Dropped {len(cnarr) - len(filtered_cn)} / {len(cnarr)} bins"
        if cnarr["chromosome"].iat[0] == cnarr["chromosome"].iat[-1]:
            msg += " on chromosome " + str(cnarr["chromosome"].iat[0])
        logging.info(msg)
    if not len(filtered_cn):
        return filtered_cn

    seg_out = ""
    match method:
        case "haar":
            # With variants, segment jointly on depth+BAF (breakpoint union) to
            # catch copy-neutral LOH; without, depth-only (byte-identical).
            segarr = haar.segment_haar(filtered_cn, threshold, variants)  # type: ignore[arg-type]
        case "none":
            segarr = none.segment_none(filtered_cn)
        case "cbs":
            # Run the R script (DNAcopy CBS) to calculate copy number segments
            rscript = cbs.CBS_RSCRIPT

            filtered_cn["start"] += 1  # Convert to 1-indexed coordinates for R
            with tempfile.NamedTemporaryFile(suffix=".cnr", mode="w+t") as tmp:
                # TODO tabio.write(filtered_cn, tmp, 'seg')
                filtered_cn.data.to_csv(
                    tmp, index=False, sep="\t", float_format="%.6g", mode="w+t"
                )
                tmp.flush()
                with core.temp_write_text(rscript, mode="w+t") as script_fname:
                    # Pass run parameters as positional command-line arguments
                    # (read via commandArgs in the R script) instead of
                    # interpolating them into the script source, so a sample ID
                    # or tempfile path containing quotes or backslashes can't
                    # produce invalid R. call_quiet runs Rscript without a shell,
                    # so the values need no escaping.
                    seg_out = core.call_quiet(
                        rscript_path,
                        "--no-restore",
                        "--no-environ",
                        script_fname,
                        tmp.name,
                        str(cnarr.sample_id),
                        str(threshold),
                        "TRUE" if smooth_cbs else "FALSE",
                    )
            # Convert R dataframe contents (SEG) to a proper CopyNumArray
            # NB: Automatically shifts 'start' back from 1- to 0-indexed
            segarr = tabio.read(StringIO(seg_out.decode()), "seg", into=CNA)  # type: ignore[arg-type,assignment]
        case "hmm" | "hmm-tumor" | "hmm-germline":
            segarr = hmm.segment_hmm(
                filtered_cn,
                method,
                diploid_parx_genome,
                int(threshold) if threshold is not None else None,
                variants,
            )
        case _:
            raise ValueError(f"Unknown method {method!r}")

    segarr.meta = cnarr.meta.copy()
    if not len(segarr):
        # No segments produced -- e.g. every bin in this arm was unsegmentable
        # in R (missing chromosome or non-finite start), so DNAcopy emitted an
        # empty SEG table. Return the empty array rather than crashing in
        # variant re-segmentation or transfer_fields, which both assume at
        # least one segment (#868).
        if save_dataframe:
            return segarr, seg_out
        return segarr
    if variants and method not in _HMM_METHODS:
        # ('hmm*' models BAF jointly in segment_hmm, so it's excluded here.)
        if method != "haar":
            # CBS and 'none': re-segment the variant allele
            # freqs within each segment using a 2-state Viterbi HMM on mirrored
            # BAF (per commit 692d5a5).
            # TODO train on all segments together
            logging.info("Re-segmenting on variant allele frequency")
            newsegs = [
                hmm.variants_in_segment(subvarr, segment)
                for segment, subvarr in variants.by_ranges(segarr)
            ]
            segarr = segarr.as_dataframe(pd.concat(newsegs))
        # else: 'haar' already incorporated BAF via its depth+BAF breakpoint
        # union, so no re-segmentation -- just attach the per-segment BAF below.
        segarr["baf"] = variants.baf_by_ranges(segarr)

    segarr = transfer_fields(segarr, cnarr)
    if save_dataframe:
        return segarr, seg_out
    return segarr


def drop_outliers(cnarr: CNA, width: int, factor: int) -> CNA:
    """Drop outlier bins with log2 ratios too far from the trend line.

    Outliers are the log2 values `factor` times the 90th quantile of absolute
    deviations from the rolling average, within a window of given `width`. The
    90th quantile is about 1.97 standard deviations if the log2 values are
    Gaussian, so this is similar to calling outliers `factor` * 1.97 standard
    deviations from the rolling mean. For a window size of 50, the breakdown
    point is 2.5 outliers within a window, which is plenty robust for our needs.
    """
    if not len(cnarr):
        return cnarr
    outlier_mask = np.concatenate(
        [
            smoothing.rolling_outlier_quantile(subarr["log2"], width, 0.95, factor)
            for _chrom, subarr in cnarr.by_chromosome()
        ]
    )
    n_outliers = outlier_mask.sum()
    if n_outliers:
        logging.info(
            "Dropped %d outlier bins:\n%s%s",
            n_outliers,
            cnarr[outlier_mask].data.head(20),
            "\n..." if n_outliers > 20 else "",
        )
    return cnarr[~outlier_mask]


def transfer_fields(
    segments: CNA, cnarr: CNA, ignore: tuple[str, ...] = params.IGNORE_GENE_NAMES
) -> CNA:
    """Map gene names, weights, depths from `cnarr` bins to `segarr` segments.

    Segment gene name is the comma-separated list of bin gene names. Segment
    weight is the sum of bin weights, and depth is the (weighted) mean of bin
    depths. Bins with NaN weights are excluded from the sums and averages.

    Also: Post-process segmentation output.

    1. Ensure every chromosome has at least one segment.
    2. Ensure first and last segment ends match 1st/last bin ends
       (but keep log2 as-is).

    """

    def make_null_segment(chrom, orig_start, orig_end):
        """Closes over 'segments'."""
        vals = {
            "chromosome": chrom,
            "start": orig_start,
            "end": orig_end,
            "gene": "-",
            "depth": 0.0,
            "log2": 0.0,
            "probes": 0.0,
            "weight": 0.0,
        }
        row_vals = tuple(vals[c] for c in segments.data.columns)
        return row_vals

    if not len(cnarr):
        # This Should Never Happen (TM)
        # raise RuntimeError("No bins for:\n" + str(segments.data))
        logging.warning("No bins for:\n%s", segments.data)
        return segments

    # Adjust segment endpoints to cover the chromosome arm's original bins
    # (Stretch first and last segment endpoints to match first/last bins)
    bins_chrom = cnarr.chromosome.iat[0]
    bins_start = cnarr.start.iat[0]
    bins_end = cnarr.end.iat[-1]
    if not len(segments):
        # All bins in this chromosome arm were dropped: make a dummy segment
        return make_null_segment(bins_chrom, bins_start, bins_end)  # type: ignore[return-value,no-any-return]
    # Avoid chained assignment by directly modifying the underlying DataFrame
    segments.data.loc[segments.data.index[0], "start"] = bins_start
    segments.data.loc[segments.data.index[-1], "end"] = bins_end

    # Aggregate segment depths, weights, gene names
    # ENH refactor so that np/CNA.data access is encapsulated in skgenome
    ignore += params.ANTITARGET_ALIASES
    assert bins_chrom == segments.chromosome.iat[0]
    cdata = cnarr.data.reset_index()
    if "depth" not in cdata.columns:
        cdata["depth"] = np.exp2(cnarr["log2"].to_numpy())
    bin_genes = cdata["gene"].to_numpy()
    bin_weights = cdata["weight"].to_numpy() if "weight" in cdata.columns else None
    bin_depths = cdata["depth"].to_numpy()
    seg_genes = ["-"] * len(segments)
    seg_weights = np.zeros(len(segments))
    seg_depths = np.zeros(len(segments))

    for i, bin_idx in enumerate(iter_slices(cdata, segments.data, "outer", False)):
        if bin_weights is not None:
            wt = bin_weights[bin_idx]
            # Use nansum so NaN weights don't propagate into .cns output
            seg_wt = float(np.nansum(wt))
            if seg_wt > 0:
                valid = ~np.isnan(wt)
                seg_dp = np.average(bin_depths[bin_idx][valid], weights=wt[valid])
            else:
                seg_dp = 0.0
        else:
            bin_count = len(cdata.iloc[bin_idx])
            seg_wt = float(bin_count)
            seg_dp = bin_depths[bin_idx].mean()
        seg_gn = join_strings(bin_genes[bin_idx], ignore=ignore)
        seg_genes[i] = seg_gn
        seg_weights[i] = seg_wt
        seg_depths[i] = seg_dp

    segments.data = segments.data.assign(
        gene=seg_genes, weight=seg_weights, depth=seg_depths
    )
    return segments
