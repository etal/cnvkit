CBS_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by CBS.
# Input: log2 coverage data in CNVkit's tabular format
# Output: the CBS data table (SEG)

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


