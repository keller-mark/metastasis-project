import pandas as pd
import numpy as np

def cellphonedb_orthologs(o_df, gi_df, pi_df, pc_df, ic_df):
  
  o_df = o_df.drop(columns=[
      "Gene stable ID version",
      "Transcript stable ID version",
      "Mouse protein or transcript stable ID"
  ])
  o_df = o_df.rename(columns={
      "Gene stable ID": "ensembl",
      "Transcript stable ID": "ensembl_transcript",
      "Mouse gene stable ID": "ensembl_gene_mouse",
      "Mouse gene name": "gene_mouse",
      "Mouse orthology confidence [0 low, 1 high]": "mouse_orthology_confidence"
  })
  #o_df = o_df.loc[o_df["mouse_orthology_confidence"] == 1]
  o_df = o_df.drop(columns=["ensembl_transcript"])
  o_df = o_df.drop_duplicates()
  
  def col_append(df, s):
      return df.rename(columns=dict(zip(df.columns.values.tolist(), [f"{col}{s}" for col in df.columns.values.tolist()])))
  
  
  gi_a_df = col_append(gi_df, "_a")
  gi_b_df = col_append(gi_df, "_b")
  
  o_a_df = col_append(o_df, "_a")
  o_b_df = col_append(o_df, "_b")
  
  ic_gi_df = ic_df.merge(gi_a_df, how="inner", left_on="partner_a", right_on="uniprot_a")
  ic_gi_df = ic_gi_df.merge(gi_b_df, how="inner", left_on="partner_b", right_on="uniprot_b")
  ic_gi_df = ic_gi_df.drop(columns=["uniprot_a", "uniprot_b"])
  
  ic_gi_o_df = ic_gi_df.merge(o_a_df, how="inner", left_on="ensembl_a", right_on="ensembl_a")
  ic_gi_o_df = ic_gi_o_df.merge(o_b_df, how="inner", left_on="ensembl_b", right_on="ensembl_b")

  return ic_gi_o_df
  
if __name__ == "__main__":
  o_df = pd.read_csv(snakemake.input['orthologs'], sep='\t', skiprows=1)
  gi_df = pd.read_csv(snakemake.input['cpdb_gene_input'])
  pi_df = pd.read_csv(snakemake.input['cpdb_protein_input'])
  pc_df = pd.read_csv(snakemake.input['cpdb_protein_curated'])
  ic_df = pd.read_csv(snakemake.input['cpdb_interaction_curated'])
  
  df = cellphonedb_orthologs(o_df, gi_df, pi_df, pc_df, ic_df)
  df.to_csv(snakemake.output[0], sep='\t', index=False)