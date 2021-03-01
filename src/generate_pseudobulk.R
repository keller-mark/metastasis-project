library(muscat)
library(Seurat)
library(SingleCellExperiment)
library(sceasy)
library(reticulate)

use_condaenv('lr-env')


# Load the Seurat object (into a variable named `tiss`)
load(snakemake@input[[1]])

tiss <- UpdateSeuratObject(tiss)
tiss_sce <- as.SingleCellExperiment(tiss)
pb <- aggregateData(tiss.sce, assay = "counts", fun = "sum", by = c("cell_ontology_id"))
pb_assays <- as.list(assays(pb))
names(pb_assays) <- c("pseudobulk")

pb_temp <- SingleCellExperiment(
  assays = pb_assays
)

assays(pb) <- assays(pb_temp)

sceasy::convertFormat(pb, from = "sce", to="anndata", outFile=snakemake@output[[1]])
