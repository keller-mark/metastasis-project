# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html#5_Estimating_parameters

library(Seurat)
library(splatter)
library(scater)
library(jsonlite)

set.seed(2445)

# Load the Seurat object (into a variable named `tiss`)
load(snakemake@input[[1]])

tissue <- snakemake@params[['tissue']]

tiss <- UpdateSeuratObject(tiss)
sparse_matrix <- slot(slot(tiss, "assays")[['RNA']], "counts")

# Each row is a gene, each column is a cell.
dense_matrix <- as.matrix(sparse_matrix)

# Estimate parameters from data
params <- simpleEstimate(dense_matrix)

params_list <- list()

for(k in slotNames(params)) {
  params_list[[k]] <- unbox(slot(params, k))
}

params_json <- toJSON(params_list)
write(params_json, snakemake@output[[1]])