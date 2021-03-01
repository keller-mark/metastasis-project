from anndata import read_h5ad, AnnData
import scanpy as sc
import numpy as np
import pandas as pd
from os.path import join
import json

def sample_bulk(adata, splatter):
  
  # The number of genes.
  n_genes = splatter["nGenes"]
  # The number of cells.
  n_cells = splatter["nCells"]
  # The shape parameter for the mean gamma distribution.
  mean_shape = splatter["mean.shape"]
  # The rate parameter for the mean gamma distribution.
  mean_rate = splatter["mean.rate"]
  # The dispersion parameter for the counts negative binomial distribution.
  count_disp = splatter["count.disp"]
  
  
  
  X = adata.X
  
  
  return adata
  

if __name__ == "__main__":
  adata = read_h5ad(snakemake.input['adata'])
  with open(snakemake.input['splatter']) as f:
    splatter = json.load(f)
  
  adata = sample_bulk(adata, splatter)
  
  adata.write(snakemake.output[0])
