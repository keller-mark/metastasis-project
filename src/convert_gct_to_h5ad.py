from cmapPy.pandasGEXpress.parse import parse as gct_parse
from anndata import AnnData
import scanpy as sc
import numpy as np
import pandas as pd

def convert_gct_to_h5ad(gct_obj):
  adata = AnnData(X=gct_obj.data_df.transpose(), obs=gct_obj.col_metadata_df, var=gct_obj.row_metadata_df)
  return adata
  
if __name__ == "__main__":
  gct_obj = gct_parse(snakemake.input[0])
  adata = convert_gct_to_h5ad(gct_obj)
  adata.write(snakemake.output[0])
  