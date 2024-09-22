CBS_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by CBS.
# Input: log2 coverage data in CNVkit's tabular format
# Output: the CBS data table (SEG)

library('DNAcopy')

write("Loading probe coverages into a data frame", stderr())
tbl = read.delim("%(probes_fname)s")
# Filter the original .cnr to valid rows (usually they're all valid)
if (!is.null(tbl$weight) && (sum(!is.na(tbl$weight)) > 0)) {
    # Drop any null or 0-weight bins
    tbl = tbl[!is.na(tbl$weight) & (tbl$weight > 0),]
} else {
    # No bin weight information; set all equal weights
    tbl$weight = 1.0
}
if ((sum(is.na(tbl$log2)) > 0)) {
    # Drop any bins with null/missing log2 ratios
    tbl = tbl[!is.na(tbl$log2),]
}

cna = CNA(cbind(tbl$log2), tbl$chromosome, tbl$start,
          data.type="logratio", presorted=T)

write("Segmenting the probe data", stderr())
set.seed(0xA5EED)

# additional smoothing (if --smooth-cbs provided)
if (%(smooth_cbs)g) {
    write("Performing smoothing of the data", stderr())
    cna = smooth.CNA(cna)
}

fit = segment(cna, weights=tbl$weight, alpha=%(threshold)g)

write("Setting segment endpoints to original bin start/end positions", stderr())
write("and recalculating segment means with bin weights", stderr())
for (idx in 1:nrow(fit$output)) {
    if (!is.na(fit$segRows$startRow[idx])) {
        start_bin = fit$segRows$startRow[idx]
        end_bin = fit$segRows$endRow[idx]
        fit$output$loc.start[idx] = tbl$start[start_bin]
        fit$output$loc.end[idx] = tbl$end[end_bin]
        fit$output$seg.mean[idx] = weighted.mean(tbl$log2[start_bin:end_bin],
                                                 tbl$weight[start_bin:end_bin])
    }
}

write("Printing the CBS table to standard output", stderr())
fit$output$ID = '%(sample_id)s'
write.table(fit$output, '', sep='\t', row.names=FALSE)
"""
