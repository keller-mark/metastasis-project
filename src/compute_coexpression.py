import pandas as pd
import numpy as np
from anndata import read_h5ad

def compute_coexpression(pb_adata, mm_counts_df, mm_samples_df, cpdb_orthologs_df):
  
  
  
if __name__ == "__main__":
  pb_adata = read_h5ad(snakemake.input['tm_pseudobulk'])
  mm_counts_df = pd.read_csv(snakemake.input['mm_counts'])
  mm_samples_df = pd.read_csv(snakemake.input['mm_samples'])
  cpdb_orthologs_df = pd.read_csv(snakemake.input['cpdb_orthologs'])
  
  
  df = compute_coexpression(pb_adata, mm_counts_df, mm_samples_df, cpdb_orthologs_df)
  df.to_csv(snakemake.output[0], sep='\t', index=False)