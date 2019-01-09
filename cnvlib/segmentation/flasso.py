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

write(paste("Segmenting", levels(tbl$chromosome)), stderr())
fit <- cghFLasso(tbl$log2, FDR=%(threshold)g, chromosome=tbl$chromosome)

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
