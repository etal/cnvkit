"""Segmentation of copy number values."""
from __future__ import absolute_import, division
import math
import os.path
import tempfile

import numpy as np
import pandas as pd

from .. import core, ngfrills, params
from ..cnary import CopyNumArray as CNA

from Bio._py3k import StringIO


def do_segmentation(probes_fname, save_dataframe, method, threshold=None,
                    rlibpath=None):
    """Infer copy number segments from the given coverage table."""
    probes = CNA.read(probes_fname)
    filtered_probes = probes.drop_low_coverage()
    if method == 'haar':
        from . import haar
        return haar.segment_haar(filtered_probes)

    # Run R scripts to calculate copy number segments
    if method == 'cbs':
        rscript = CBS_RSCRIPT
        threshold = threshold or 0.0001
    elif method == 'flasso':
        rscript = FLASSO_RSCRIPT
        threshold = threshold or 0.005
    else:
        raise ValueError("Unknown method %r" % method)

    with tempfile.NamedTemporaryFile(suffix='.cnr') as tmp:
        filtered_probes.data.to_csv(tmp, index=False, sep='\t',
                                    float_format='%.6g')
        tmp.flush()
        script_strings = {
            'probes_fname': tmp.name,
            'sample_id': probes.sample_id,
            'threshold': threshold,
            'rlibpath': ('.libPaths(c("%s"))' % rlibpath if rlibpath else ''),
        }
        with ngfrills.temp_write_text(rscript % script_strings) as script_fname:
            seg_out = ngfrills.call_quiet('Rscript', script_fname)
        # ENH: run each chromosome separately
        # ENH: run each chrom. arm separately (via knownsegs)
    seg_pset = probes.as_dataframe(seg2cns(seg_out))

    if method == 'flasso':
        seg_pset = squash_segments(seg_pset)

    seg_pset = repair_segments(seg_pset, probes)

    if save_dataframe:
        return seg_pset, seg_out
    else:
        return seg_pset


def seg2cns(seg_text):
    """Convert R dataframe contents (SEG) to our native tabular format.

    Return a pandas.Dataframe with CNA columns.
    """
    table = pd.read_table(StringIO(seg_text), comment='[')
    if len(table.columns) == 6:
        table.columns = ["sample_id", "chromosome", "start", "end", "probes",
                         "log2"]
    elif len(table.columns) == 5:
        table.columns = ["sample_id", "chromosome", "start", "end", "log2"]
    else:
        raise ValueError("Segmentation output is not valid SEG format:\n"
                        + seg_text)
    del table["sample_id"]
    table["start"] = [int(math.ceil(float(val))) for val in table["start"]]
    table["end"] = [int(math.ceil(float(val))) for val in table["end"]]
    table["gene"] = '-'
    return table


def squash_segments(seg_pset):
    """Combine contiguous segments."""
    curr_chrom = None
    curr_start = None
    curr_end = None
    curr_val = None
    curr_cnt = 0
    squashed_rows = []
    for row in seg_pset:
        if row['chromosome'] == curr_chrom and row['log2'] == curr_val:
            # Continue the current segment
            curr_end = row['end']
            curr_cnt += 1
        else:
            # Segment break
            # Finish the current segment
            if curr_cnt:
                squashed_rows.append((curr_chrom, curr_start, curr_end,
                                      ('G' if curr_val >= 0. else 'L'),
                                      curr_val, curr_cnt))
            # Start a new segment
            curr_chrom = row['chromosome']
            curr_start = row['start']
            curr_end = row['end']
            curr_val = row['log2']
            curr_cnt = 1
    # Remainder
    squashed_rows.append((curr_chrom, curr_start, curr_end,
                          ('G' if curr_val >= 0. else 'L'),
                          curr_val, curr_cnt))
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
            null_segment = (chrom, orig_start, orig_end, "_", 0.0, 0)
            extra_segments.append(null_segment)
    if extra_segments:
        segments.concat(segments.as_rows(extra_segments))
    # Adjust segment values
    for i, (_seg, subprobes) in enumerate(orig_probes.by_segment(segments)):
        # Segment name is the comma-separated list of bin gene names
        subgenes = [g for g in pd.unique(subprobes['gene'])
                    if g not in ('Background', 'CGH', '-')]
        segments[i, 'gene'] = ",".join(subgenes) if subgenes else '-'
        # ENH: Recalculate segment means here instead of in R
    segments.sort_columns()
    return segments


CBS_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by CBS.
# Input: log2 coverage data in Nexus 'basic' format
# Output: the CBS data table

%(rlibpath)s
library('PSCBS') # Requires: R.utils, R.oo, R.methodsS3

write("Loading probe coverages into a data frame", stderr())
tbl = read.delim("%(probes_fname)s")
chrom_rle = rle(as.character(tbl$chromosome))
chrom_names = chrom_rle$value
chrom_lengths = chrom_rle$lengths
chrom_ids = rep(1:length(chrom_names), chrom_lengths)
if (is.null(tbl$weight)) {
    cna = data.frame(chromosome=chrom_ids, x=tbl$start, y=tbl$log2)
} else {
    cna = data.frame(chromosome=chrom_ids, x=tbl$start, y=tbl$log2, w=tbl$weight)
}

write("Pre-processing the probe data for segmentation", stderr())
# Find and exclude the centromere of each chromosome
largegaps = findLargeGaps(cna, minLength=1e6)
if (is.null(largegaps)) {
    knownsegs = NULL
} else {
    # Choose the largest gap in each chromosome and only omit that
    rows_to_keep = c()
    for (i in 1:length(chrom_names)) {
        curr_chrom_mask = (largegaps$chromosome == i)
        if (sum(curr_chrom_mask)) {
            best = which(
                curr_chrom_mask &
                (largegaps$length == max(largegaps[curr_chrom_mask,]$length))
            )
            rows_to_keep = c(rows_to_keep, best)
        }
    }
    knownsegs = gapsToSegments(largegaps[rows_to_keep,])
}

write("Segmenting the probe data", stderr())
fit = segmentByCBS(cna, alpha=%(threshold)g, undo=0, min.width=2,
                   joinSegments=FALSE, knownSegments=knownsegs, seed=0xA5EED)

write("Setting segment endpoints to original bin start/end positions", stderr())
write("and recalculating segment means with bin weights", stderr())
for (idx in 1:nrow(fit$output)) {
    if (!is.na(fit$segRows$startRow[idx])) {
        start_bin = fit$segRows$startRow[idx]
        end_bin = fit$segRows$endRow[idx]
        fit$output$start[idx] = tbl$start[start_bin]
        fit$output$end[idx] = tbl$end[end_bin]
        fit$output$mean[idx] = weighted.mean(tbl$log2[start_bin:end_bin],
                                             tbl$weight[start_bin:end_bin])
    }
}

write("Restoring the original chromosome names", stderr())
fit$output$sampleName = '%(sample_id)s'
out = na.omit(fit$output) # Copy for lookup in the loop
out2 = na.omit(fit$output) # Copy to modify
for (i in 1:length(chrom_names)) {
    out2[out$chromosome == i,]$chromosome = chrom_names[i]
}

write("Printing the CBS table to standard output", stderr())
write.table(out2, '', sep='\t', row.names=FALSE)
"""


FLASSO_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by CBS.
# Input: log2 coverage data in Nexus 'basic' format
# Output: the CBS data table

%(rlibpath)s
library('cghFLasso')

tbl <- read.delim("%(probes_fname)s")
# Ignore low-coverage probes
positions <- (tbl$start + tbl$end) * 0.5

write("Segmenting the probe data", stderr())
fit <- cghFLasso(tbl$log2, FDR=%(threshold)g)

# Reformat the output table as SEG
outtable <- data.frame(sample="%(sample_id)s",
                       chromosome=tbl$chromosome,
                       start=tbl$start,
                       end=tbl$end,
                       nprobes=1,
                       value=fit$Esti.CopyN)

write("Printing the segment table to standard output", stderr())
write.table(outtable, '', sep='\t', row.names=FALSE)
"""
