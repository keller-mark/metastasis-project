import pandas as pd
import numpy as np

def cellphonedb_orthologs(o_df, gi_df, pi_df, pc_df):
  
  print(o_df.head())
  return df
  
if __name__ == "__main__":
  o_df = pd.read_csv(snakemake.input['orthologs'], sep='\t', skiprows=1)
  gi_df = pd.read_csv(snakemake.input['cpdb_gene_input'])
  pi_df = pd.read_csv(snakemake.input['cpdb_prot_input'])
  pc_df = pd.read_csv(snakemake.input['cpdb_prot_curated'])
  
  df = cellphonedb_orthologs(o_df, gi_df, pi_df, pc_df)
  df.to_csv(snakemake.output[0], sep='\t', index=False)