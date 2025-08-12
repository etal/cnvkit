"""Segmentation of copy number values."""

from __future__ import annotations
import locale
import logging
import tempfile
from io import StringIO
from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd
from skgenome import tabio
from skgenome.intersect import iter_slices

from .. import core, parallel, params, smoothing
from ..cnary import CopyNumArray as CNA
from ..segfilters import squash_by_groups
from . import cbs, flasso, haar, none

if TYPE_CHECKING:
    from cnvlib.vary import VariantArray

try:
    from . import hmm

    HMM_METHODS = ("hmm", "hmm-tumor", "hmm-germline")
except ImportError as e:
    hmm = None
    HMM_METHODS = ()
    HMM_IMPORT_ERROR = str(e)

SEGMENT_METHODS = ("cbs", "flasso", "haar", "none", *HMM_METHODS)


def do_segmentation(
    cnarr: CNA,
    method: str,
    diploid_parx_genome: None = None,
    threshold: Optional[float] = None,
    variants: Optional[VariantArray] = None,
    skip_low: bool = False,
    skip_outliers: int = 10,
    min_weight: int = 0,
    save_dataframe: bool = False,
    rscript_path: str = "Rscript",
    processes: int = 1,
    smooth_cbs: bool = False,
) -> CNA:
    """Infer copy number segments from the given coverage table."""
    if method not in SEGMENT_METHODS:
        if method in ("hmm", "hmm-tumor", "hmm-germline") and not HMM_METHODS:
            raise ImportError(
                f"HMM segmentation method '{method}' requires pomegranate >= 1.0.0. "
                f"Install with: pip install cnvkit[hmm] or pip install pomegranate>=1.0.0. "
                f"Error: {HMM_IMPORT_ERROR}"
            )
        raise ValueError(
            "'method' must be one of: "
            + ", ".join(SEGMENT_METHODS)
            + "; got: "
            + repr(method)
        )

    if not threshold:
        threshold = {
            "cbs": 0.0001,
            "flasso": 0.0001,
            "haar": 0.0001,
        }.get(method)
    msg = "Segmenting with method " + repr(method)
    if threshold is not None:
        if method.startswith("hmm"):
            msg += f", smoothing window size {threshold},"
        else:
            msg += f", significance threshold {threshold},"
    msg += f" in {processes} processes"
    logging.info(msg)

    # NB: parallel cghFLasso segfaults in R ('memory not mapped'),
    # even when run on a single chromosome
    if method == "flasso" or method.startswith("hmm"):
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
        )
        if save_dataframe:
            cna, rstr = cna
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
            rets, r_dframe_strings = zip(*rets)
            # Strip the header line from all but the first dataframe, then combine
            r_dframe_strings = map(_to_str, r_dframe_strings)
            rstr = [next(r_dframe_strings)]
            rstr.extend(r[r.index("\n") + 1 :] for r in r_dframe_strings)
            rstr = "".join(rstr)
        cna = cnarr.concat(rets)

    cna.sort_columns()
    if save_dataframe:
        return cna, rstr
    return cna


def _to_str(s, enc=locale.getpreferredencoding()):  # noqa: B008
    # Evaluate encoding once at function definition for performance
    if isinstance(s, bytes):
        return s.decode(enc)
    return s


def _ds(
    args: tuple[CNA, str, None, float, None, bool, int, int, bool, str, bool],
) -> CNA:
    """Wrapper for parallel map"""
    return _do_segmentation(*args)


def _do_segmentation(
    cnarr: CNA,
    method: str,
    diploid_parx_genome: None,
    threshold: Optional[float],
    variants: Optional[VariantArray] = None,
    skip_low: bool = False,
    skip_outliers: int = 10,
    min_weight: int = 0,
    save_dataframe: bool = False,
    rscript_path: str = "Rscript",
    smooth_cbs: bool = False,
) -> CNA:
    """Infer copy number segments from the given coverage table."""
    if not len(cnarr):
        return cnarr

    filtered_cn = cnarr.copy()
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
    if method == "haar":
        segarr = haar.segment_haar(filtered_cn, threshold)

    elif method == "none":
        segarr = none.segment_none(filtered_cn)

    elif method.startswith("hmm"):
        segarr = hmm.segment_hmm(
            filtered_cn, method, diploid_parx_genome, threshold, variants
        )

    elif method in ("cbs", "flasso"):
        # Run R scripts to calculate copy number segments
        rscript = {
            "cbs": cbs.CBS_RSCRIPT,
            "flasso": flasso.FLASSO_RSCRIPT,
        }[method]

        filtered_cn["start"] += 1  # Convert to 1-indexed coordinates for R
        with tempfile.NamedTemporaryFile(suffix=".cnr", mode="w+t") as tmp:
            # TODO tabio.write(filtered_cn, tmp, 'seg')
            filtered_cn.data.to_csv(
                tmp, index=False, sep="\t", float_format="%.6g", mode="w+t"
            )
            tmp.flush()
            script_strings = {
                "probes_fname": tmp.name,
                "sample_id": cnarr.sample_id,
                "threshold": threshold,
                "smooth_cbs": smooth_cbs,
            }
            with core.temp_write_text(
                rscript % script_strings, mode="w+t"
            ) as script_fname:
                seg_out = core.call_quiet(
                    rscript_path, "--no-restore", "--no-environ", script_fname
                )
        # Convert R dataframe contents (SEG) to a proper CopyNumArray
        # NB: Automatically shifts 'start' back from 1- to 0-indexed
        segarr = tabio.read(StringIO(seg_out.decode()), "seg", into=CNA)
        if method == "flasso":
            # Merge adjacent bins with same log2 value into segments
            if "weight" in filtered_cn:
                segarr["weight"] = filtered_cn["weight"]
            else:
                segarr["weight"] = 1.0
            segarr = squash_by_groups(segarr, segarr["log2"], by_arm=True)

    else:
        raise ValueError(f"Unknown method {method!r}")

    segarr.meta = cnarr.meta.copy()
    if variants and not method.startswith("hmm"):
        # Re-segment the variant allele freqs within each segment
        # TODO train on all segments together
        logging.info("Re-segmenting on variant allele frequency")
        newsegs = [
            hmm.variants_in_segment(subvarr, segment)
            for segment, subvarr in variants.by_ranges(segarr)
        ]
        segarr = segarr.as_dataframe(pd.concat(newsegs))
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
    segments: CNA, cnarr: CNA, ignore: tuple[str, str, str] = params.IGNORE_GENE_NAMES
) -> CNA:
    """Map gene names, weights, depths from `cnarr` bins to `segarr` segments.

    Segment gene name is the comma-separated list of bin gene names. Segment
    weight is the sum of bin weights, and depth is the (weighted) mean of bin
    depths.

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
        return make_null_segment(bins_chrom, bins_start, bins_end)
    segments.start.iat[0] = bins_start
    segments.end.iat[-1] = bins_end

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
            seg_wt = bin_weights[bin_idx].sum()
            if seg_wt > 0:
                seg_dp = np.average(bin_depths[bin_idx], weights=bin_weights[bin_idx])
            else:
                seg_dp = 0.0
        else:
            bin_count = len(cdata.iloc[bin_idx])
            seg_wt = float(bin_count)
            seg_dp = bin_depths[bin_idx].mean()
        subgenes = [g for g in pd.unique(bin_genes[bin_idx]) if g not in ignore]
        seg_gn = ",".join(subgenes) if subgenes else "-"
        seg_genes[i] = seg_gn
        seg_weights[i] = seg_wt
        seg_depths[i] = seg_dp

    segments.data = segments.data.assign(
        gene=seg_genes, weight=seg_weights, depth=seg_depths
    )
    return segments
