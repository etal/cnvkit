"""Segmentation of copy number values."""
from __future__ import absolute_import, division, print_function
import logging
import math
import os.path
import tempfile

import numpy as np
import pandas as pd

from .. import core, params, smoothing, tabio, vary
from ..cnary import CopyNumArray as CNA
from . import cbs, flasso, haar

from Bio._py3k import StringIO


def do_segmentation(cnarr, method, threshold=None, variants=None,
                    skip_low=False, skip_outliers=10,
                    save_dataframe=False, rlibpath=None):
    """Infer copy number segments from the given coverage table."""
    filtered_cn = cnarr
    if skip_low:
        before = len(filtered_cn)
        filtered_cn = filtered_cn.drop_low_coverage()
        logging.info("Dropped %d low-coverage bins", before - len(filtered_cn))
    if skip_outliers:
        filtered_cn = drop_outliers(filtered_cn, 50, skip_outliers)

    if method == 'haar':
        threshold = threshold or 0.001
        segarr = haar.segment_haar(filtered_cn, threshold)
        segarr['gene'], segarr['weight'] = transfer_names_weights(segarr, cnarr)

    elif method in ('cbs', 'flasso'):
        # Run R scripts to calculate copy number segments
        if method == 'cbs':
            rscript = cbs.CBS_RSCRIPT
            threshold = threshold or 0.0001
        elif method == 'flasso':
            rscript = flasso.FLASSO_RSCRIPT
            threshold = threshold or 0.005

        with tempfile.NamedTemporaryFile(suffix='.cnr', mode="w+t") as tmp:
            filtered_cn.data.to_csv(tmp, index=False, sep='\t',
                                    float_format='%.6g', mode="w+t")
            tmp.flush()
            script_strings = {
                'probes_fname': tmp.name,
                'sample_id': cnarr.sample_id,
                'threshold': threshold,
                'rlibpath': ('.libPaths(c("%s"))' % rlibpath if rlibpath else ''),
            }
            with core.temp_write_text(rscript % script_strings,
                                          mode="w+t") as script_fname:
                seg_out = core.call_quiet('Rscript', '--vanilla', script_fname)
            # ENH: run each chromosome separately
            # ENH: run each chrom. arm separately (via knownsegs)
        # Convert R dataframe contents (SEG) to a proper CopyNumArray
        segarr = tabio.read(StringIO(seg_out.decode()), "seg",
                            sample_id=cnarr.sample_id)
        if method == 'flasso':
            segarr = squash_segments(segarr)
        segarr = repair_segments(segarr, cnarr)

    else:
        raise ValueError("Unknown method %r" % method)

    if variants:
        # Re-segment the variant allele freqs within each segment
        newsegs = [haar.variants_in_segment(subvarr, segment, 0.01 * threshold)
                   for segment, subvarr in variants.by_ranges(segarr)]
        segarr = segarr.as_dataframe(pd.concat(newsegs))
        segarr.sort_columns()
        # TODO fix ploidy on allosomes
        allelics = vary._allele_specific_copy_numbers(segarr, variants)
        segarr.data = pd.concat([segarr.data, allelics], axis=1, copy=False)

    segarr['gene'], segarr['weight'] = transfer_names_weights(segarr, cnarr)

    if save_dataframe:
        return segarr, seg_out
    else:
        return segarr


def drop_outliers(cnarr, width, factor):
    """Drop outlier bins with log2 ratios too far from the trend line.

    Outliers are the log2 values `factor` times the 90th quantile of absolute
    deviations from the rolling average, within a window of given `width`. The
    90th quantile is about 1.97 standard deviations if the log2 values are
    Gaussian, so this is similar to calling outliers `factor` * 1.97 standard
    deviations from the rolling mean. For a window size of 50, the breakdown
    point is 2.5 outliers within a window, which is plenty robust for our needs.
    """
    outlier_mask = np.concatenate([
        smoothing.rolling_outlier_quantile(subarr['log2'], width, .95, factor)
        for _chrom, subarr in cnarr.by_chromosome()])
    n_outliers = outlier_mask.sum()
    if n_outliers:
        logging.info("Dropped %d outlier bins:\n%s%s",
                     n_outliers,
                     cnarr[outlier_mask].data.head(20),
                     "\n..." if n_outliers > 20 else "")
    else:
        logging.info("No outlier bins")
    return cnarr[~outlier_mask]


def transfer_names_weights(segments, cnarr, ignore=params.IGNORE_GENE_NAMES):
    """Copy gene names from `cnarr` to the segmented `segarr`.

    Segment name is the comma-separated list of bin gene names.
    """
    ignore += ("Background",)
    segnames = ['-'] * len(segments)
    segweights = np.zeros(len(segments))
    for i, (_seg, subprobes) in enumerate(cnarr.by_ranges(segments)):
        segweights[i] = (subprobes['weight'].sum()
                         if "weight" in cnarr
                         else len(subprobes))
        subgenes = [g for g in pd.unique(subprobes['gene']) if g not in ignore]
        if subgenes:
            segnames[i] = ",".join(subgenes)
    return segnames, segweights


def squash_segments(seg_pset):
    """Combine contiguous segments."""
    curr_chrom = None
    curr_start = None
    curr_end = None
    curr_genes = []
    curr_val = None
    curr_cnt = 0
    squashed_rows = []
    for row in seg_pset:
        if row.chromosome == curr_chrom and row.log2 == curr_val:
            # Continue the current segment
            curr_end = row.end
            curr_genes.append(row.gene)
            curr_cnt += 1
        else:
            # Segment break
            # Finish the current segment
            if curr_cnt:
                squashed_rows.append((curr_chrom, curr_start, curr_end,
                                      ",".join(pd.unique(curr_genes)),
                                      curr_val, 2 ** curr_val, curr_cnt))
            # Start a new segment
            curr_chrom = row.chromosome
            curr_start = row.start
            curr_end = row.end
            curr_genes = []
            curr_val = row.log2
            curr_cnt = 1
    # Remainder
    squashed_rows.append((curr_chrom, curr_start, curr_end,
                          ",".join(pd.unique(curr_genes)),
                          curr_val, 2** curr_val, curr_cnt))
    return seg_pset.as_rows(squashed_rows)


def repair_segments(segments, orig_probes):
    """Post-process segmentation output.

    1. Ensure every chromosome has at least one segment.
    2. Ensure first and last segment ends match 1st/last bin ends
       (but keep log2 as-is).
    3. Store probe-level gene names, comma-separated, as the segment name.
    """
    segments = segments.copy()
    extra_segments = []
    # Adjust segment endpoints on each chromosome
    for chrom, subprobes in orig_probes.by_chromosome():
        chr_seg_idx = np.where(segments.chromosome == chrom)[0]
        orig_start = subprobes[0, 'start']
        orig_end =  subprobes[len(subprobes)-1, 'end']
        if len(chr_seg_idx):
            segments[chr_seg_idx[0], 'start'] = orig_start
            segments[chr_seg_idx[-1], 'end'] = orig_end
        else:
            null_segment = (chrom, orig_start, orig_end, "-", 0.0, 0)
            extra_segments.append(null_segment)
    if extra_segments:
        segments.add(segments.as_rows(extra_segments))
    # ENH: Recalculate segment means here instead of in R
    return segments
