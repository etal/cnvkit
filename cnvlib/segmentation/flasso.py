FLASSO_RSCRIPT = """\
#!/usr/bin/env Rscript

# Calculate copy number segmentation by Fused Lasso.
# Input: log2 coverage data in CNVkit's tabular format
# Output: the CBS-style SEG data table

library('cghFLasso')

# Read parameters as positional command-line arguments rather than
# interpolating them into the script source, so sample IDs or file paths
# containing quotes or backslashes cannot break the R parser.
args = commandArgs(trailingOnly=TRUE)
probes_fname = args[1]
sample_id = args[2]
threshold = as.numeric(args[3])

tbl <- read.delim(probes_fname)

write(paste("Segmenting", levels(tbl$chromosome)), stderr())
fit <- cghFLasso(tbl$log2, FDR=threshold, chromosome=tbl$chromosome)

# Reformat the output table as SEG
outtable <- data.frame(sample=sample_id,
                       chromosome=tbl$chromosome,
                       start=tbl$start,
                       end=tbl$end,
                       nprobes=1,
                       value=fit$Esti.CopyN)

write("Printing the segment table to standard output", stderr())
write.table(outtable, '', sep='\t', row.names=FALSE)
"""
