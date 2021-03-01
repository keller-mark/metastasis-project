from biomart import BiomartServer
import pandas as pd
import numpy as np

def download_orthologs():
  # Reference: https://github.com/sebriois/biomart
  server = BiomartServer( "http://uswest.ensembl.org/biomart")
  dataset = server.datasets['hsapiens_gene_ensembl']
  
  # Reference: https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
  response = dataset.search({
    'filters': {
      'with_mmusculus_homolog': 'only'
    },
    'attributes': [
      "ensembl_gene_id",
      "ensembl_gene_id_version",
    	"ensembl_transcript_id",
    	"ensembl_transcript_id_version",
    	"mmusculus_homolog_ensembl_gene",
    	"mmusculus_homolog_associated_gene_name",
    	"mmusculus_homolog_ensembl_peptide",
    	"mmusculus_homolog_orthology_confidence",
    ]
  }, header=1)
  
  df = pd.DataFrame([ line.decode("utf-8").split("\t") for line in response.iter_lines() ])
  
  return df
  
if __name__ == "__main__":
  df = download_orthologs()
  df.to_csv(snakemake.output[0], sep='\t', index=False)