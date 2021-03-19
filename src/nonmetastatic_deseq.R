library(DESeq2)

counts_mtx <- t(as.matrix(read.csv(snakemake@input[['counts']], row.names = 1)))
conditions_df <- read.csv(snakemake@input[['conditions']], row.names = 1)

# conditions_df should have only one column
# ("condition", binary values for metastatic vs. non-metastatic) besides its index (cell line names)

dds <- DESeqDataSetFromMatrix(countData = counts_mtx, colData = conditions_df, design = ~ metastatic)

# Reference: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res),  file=snakemake@output[[1]])