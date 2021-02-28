from anndata import AnnData
import scanpy as sc
import numpy as np
import pandas as pd
from os.path import join


  
def convert_to_h5ad(X_df, ann_df, tissue):
  obs_df = ann_df.loc[ann_df["tissue"] == tissue].set_index("cell")
    
  dummy_join_cols = X_df.columns.values.tolist()[0:1]
  joined_obs_df = X_df[dummy_join_cols].join(obs_df).drop(columns=dummy_join_cols)
  
  tsne = joined_obs_df[['tissue_tSNE_1', 'tissue_tSNE_2']]
  
  drop_cols = ["Neurog3>0_raw", "Neurog3>0_scaled", 'tissue_tSNE_1', 'tissue_tSNE_2']
  joined_obs_df = joined_obs_df.drop(columns=drop_cols)
  
  adata = AnnData(X=X_df, obs=joined_obs_df, obsm={ "X_tsne": tsne.values })
  
  return adata
  

if __name__ == "__main__":
  
  X_df = pd.read_csv(snakemake.input['counts'], index_col=0).transpose()
  ann_df = pd.read_csv(snakemake.input['annots'])
  tissue = snakemake.wildcards['tissue']
  
  adata = convert_to_h5ad(X_df, ann_df, tissue)
  
  adata.write(snakemake.output[0])
  
