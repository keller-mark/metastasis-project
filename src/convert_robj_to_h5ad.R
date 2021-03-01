library(Seurat)
library(sceasy)
library(reticulate)

use_condaenv('lr-env')

# Load the Seurat object (into a variable named `tiss`)
load(snakemake@input[[1]])

tiss <- UpdateSeuratObject(tiss)

sceasy::convertFormat(tiss, from = "seurat", to="anndata", outFile=snakemake@output[[1]])
