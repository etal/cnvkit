CBS_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by CBS.
# Input: log2 coverage data in CNVkit's tabular format
# Output: the CBS data table (SEG)

%(rlibpath)s
library('DNAcopy')

write("Loading probe coverages into a data frame", stderr())
tbl = read.delim("%(probes_fname)s")
if (!is.null(tbl$weight)) {
    # Drop any 0-weight bins
    tbl = tbl[tbl$weight > 0,]
}
cna = CNA(cbind(tbl$log2), tbl$chromosome, tbl$start,
          data.type="logratio", presorted=T)

write("Segmenting the probe data", stderr())
set.seed(0xA5EED)
if (is.null(tbl$weight)) {
    fit = segment(cna, alpha=%(threshold)g)
} else {
    fit = segment(cna, weights=tbl$weight, alpha=%(threshold)g)
}

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
