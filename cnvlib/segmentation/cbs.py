CBS_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by CBS.
# Input: log2 coverage data in CNVkit's tabular format
# Output: the CBS data table (SEG)

library('DNAcopy')

# Read parameters as positional command-line arguments rather than
# interpolating them into the script source, so sample IDs or file paths
# containing quotes or backslashes cannot break the R parser.
args = commandArgs(trailingOnly=TRUE)
probes_fname = args[1]
sample_id = args[2]
threshold = as.numeric(args[3])
smooth_cbs = as.logical(args[4])

write("Loading probe coverages into a data frame", stderr())
tbl = read.delim(probes_fname)
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
# Drop bins CNA() would silently remove (missing chromosome or non-finite
# start) so the 'weights' vector below stays aligned with the probes (#868)
tbl = tbl[!is.na(tbl$chromosome) & is.finite(tbl$start),]
if (nrow(tbl) == 0) {
    # Every bin was unsegmentable; emit an empty SEG table so the caller
    # substitutes a placeholder segment instead of crashing (#868)
    empty = data.frame(ID=character(), chrom=character(), loc.start=numeric(),
                       loc.end=numeric(), num.mark=integer(), seg.mean=numeric())
    write.table(empty, '', sep='\t', row.names=FALSE)
    quit(save="no", status=0)
}

cna = CNA(cbind(tbl$log2), tbl$chromosome, tbl$start,
          data.type="logratio", presorted=T)

write("Segmenting the probe data", stderr())
set.seed(0xA5EED)

# additional smoothing (if --smooth-cbs provided)
if (smooth_cbs) {
    write("Performing smoothing of the data", stderr())
    cna = smooth.CNA(cna)
}

fit = segment(cna, weights=tbl$weight, alpha=threshold)

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
fit$output$ID = sample_id
write.table(fit$output, '', sep='\t', row.names=FALSE)
"""
