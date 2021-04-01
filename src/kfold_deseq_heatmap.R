library(viridis)
library(pheatmap)

counts_df <- read.csv(snakemake@input[['counts']], row.names = 1)
conditions_df <- read.csv(snakemake@input[['conditions']], row.names = 1)
deseq_df <- read.csv(snakemake@input[['deseq']], row.names = 1)

significant <- which(deseq_df[["padj"]] <= 0.05 & abs(deseq_df[["log2FoldChange"]]) >= 2.0)

counts_sig_df <- counts_df[,significant]
deseq_sig_df <- deseq_df[significant,c("padj", "log2FoldChange")]

conditions_df <- conditions_df[order(conditions_df$penetrance),]
deseq_sig_df <- deseq_sig_df[order(deseq_sig_df$log2FoldChange),]
counts_sig_df <- counts_sig_df[rownames(conditions_df), rownames(deseq_sig_df)]

plot <- pheatmap(
  mat = log(as.matrix(counts_sig_df)+1),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = inferno(10),
  border_color = NA,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_row = conditions_df,
  annotation_col = deseq_sig_df,
  main = paste(
    "Cell lines by differentially expressed genes",
    snakemake@wildcards[['tissue']], "fold", snakemake@wildcards[['fold']]
  )
)

save_pheatmap_pdf <- function(x, filename, width=11, height=8.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plot, snakemake@output[['heatmap_plot']])