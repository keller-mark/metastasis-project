from os.path import join
from urllib.parse import quote_plus

configfile: "config.yml"

assert("params" in config.keys())

print(config["params"])

PARAMS = config["params"]

DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
PROCESSED_DIR = join(DATA_DIR, "processed")


# Tabula Muris URLs: https://figshare.com/articles/dataset/Robject_files_for_tissues_processed_by_Seurat/5821263
TM_FACS_SEURAT = {
 'Bladder': 'https://ndownloader.figshare.com/files/13091141',
 'Brain_Non-Myeloid': 'https://ndownloader.figshare.com/files/13091513',
 'Fat': 'https://ndownloader.figshare.com/files/13091744',
 'Skin': 'https://ndownloader.figshare.com/files/13092389',
 'Trachea': 'https://ndownloader.figshare.com/files/13092410',
 'Lung': 'https://ndownloader.figshare.com/files/13092194',
 'Kidney': 'https://ndownloader.figshare.com/files/13091963',
 'Brain_Myeloid': 'https://ndownloader.figshare.com/files/13091333',
 'Pancreas': 'https://ndownloader.figshare.com/files/13092386',
 'Tongue': 'https://ndownloader.figshare.com/files/13092401',
 'Limb_Muscle': 'https://ndownloader.figshare.com/files/13092152',
 'Thymus': 'https://ndownloader.figshare.com/files/13092398',
 'Mammary_Gland': 'https://ndownloader.figshare.com/files/13092197',
 'Diaphragm': 'https://ndownloader.figshare.com/files/13091525',
 'Marrow': 'https://ndownloader.figshare.com/files/13092380',
 'Spleen': 'https://ndownloader.figshare.com/files/13092395',
 'Large_Intestine': 'https://ndownloader.figshare.com/files/13092143',
 'Liver': 'https://ndownloader.figshare.com/files/13092155',
 'Heart': 'https://ndownloader.figshare.com/files/13091957',
 'Aorta': 'https://ndownloader.figshare.com/files/13091138'
}
TM_FACS_TISSUES = list(TM_FACS_SEURAT.keys())

TM_TO_METMAP = {
  'Brain_Non-Myeloid': 'brain',
  'Lung': 'lung',
  'Liver': 'liver',
  'Marrow': 'bone',
  'Kidney': 'kidney'
}
METMAP_TISSUES = list(TM_TO_METMAP.keys())

# CellPhoneDB URLs: https://www.cellphonedb.org/downloads
CELLPHONEDB_GENE_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv"
CELLPHONEDB_PROTEIN_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/protein_input.csv"
CELLPHONEDB_COMPLEX_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/complex_input.csv"
CELLPHONEDB_INTERACTION_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/interaction_input.csv"

# MetMap URLs: https://depmap.org/metmap/data/index.html
METMAP_500_URL = "https://ndownloader.figshare.com/files/24009293"
METMAP_EXP_COUNTS_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE148283&format=file&file=GSE148283%5Fall%2Ecount%2Ecsv%2Egz"
METMAP_EXP_SAMPLES_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE148283&format=file&file=GSE148283%5Fall%2Esample%2Ecsv%2Egz"

# CCLE URLS: https://portals.broadinstitute.org/ccle/data
CCLE_EXP_URL = "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz"

# Cell and Gene Ontology URLs: http://www.obofoundry.org/ontology/cl.html
CL_OBO_URL = "http://purl.obolibrary.org/obo/cl.obo"
GO_OBO_URL = "http://purl.obolibrary.org/obo/go.obo"

# Abbreviations
# - iea: interaction_expression_aggregation (how to combine expression for partners A and B of an interaction)
# - cea: complex_expression_aggregation (how to aggregate expression for complexes which contain multiple gene products)

# Rules
rule all:
  input:
    expand(
      join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.model.h5ad"),
      tissue=METMAP_TISSUES,
      model=PARAMS["models"],
      expression_scale=PARAMS["expression_scales"],
      interaction_source=PARAMS["interaction_sources"],
      cea=PARAMS["complex_expression_aggregation"],
      iea=PARAMS["interaction_expression_aggregation"],
    ),
    expand(
      join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.mse_plot.pdf"),
      tissue=METMAP_TISSUES,
      model=PARAMS["models"],
      expression_scale=PARAMS["expression_scales"],
      interaction_source=PARAMS["interaction_sources"],
      cea=PARAMS["complex_expression_aggregation"],
      iea=PARAMS["interaction_expression_aggregation"],
    ),
    join(PROCESSED_DIR, "plots", "ensembl_orthologs.pdf"),
    join(PROCESSED_DIR, "plots", "cellphonedb_interactions.pdf"),
    join(PROCESSED_DIR, "metmap", "pca.pdf")


# Use co-expression values to build a model of metastasis potential
# for each MetMap target site.
rule build_model:
  input:
    join(PROCESSED_DIR, "coexpression", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression.h5ad"),
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue])
  output:
    model=join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.model.h5ad"),
    mse_plot=join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.mse_plot.pdf"),
    prediction_plot=join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.prediction_plot.pdf")
  notebook:
    join("src", "build_model.py.ipynb")

# Compute co-expression of CellPhoneDB interaction genes between
# each cancer cell line used in MetMap (expression from CCLE)
# and each healthy organ used in MetMap (expression from Tabula Muris)
# where the genes are linked through human-mouse orthologs from Ensembl.
rule compute_coexpression:
  input:
    tm_pseudobulk=join(RAW_DIR, "tm", "anndata", "{tissue}.pseudobulk.h5ad"),
    cl_obo=join(RAW_DIR, "ontologies", "cl.obo"),
    mm_potential=join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx"),
    ccle_exp=join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.h5ad"),
    interactions=join(PROCESSED_DIR, "cellphonedb", "gene_orthologs.tsv")
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue])
  output:
    join(PROCESSED_DIR, "coexpression", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression.h5ad")
  notebook:
    join("src", "compute_coexpression.py.ipynb")

rule cellphonedb_orthologs:
  input:
    orthologs=join(RAW_DIR, "ensembl", "human_mouse_orthologs.tsv"),
    cpdb_gene_input=join(RAW_DIR, "cellphonedb", "gene_input.csv"),
    cpdb_protein_input=join(RAW_DIR, "cellphonedb", "protein_input.csv"),
    cpdb_complex_input=join(RAW_DIR, "cellphonedb", "complex_input.csv"),
    cpdb_interaction_input=join(RAW_DIR, "cellphonedb", "interaction_input.csv")
  output:
    table=join(PROCESSED_DIR, "cellphonedb", "gene_orthologs.tsv"),
    ortholog_plot=join(PROCESSED_DIR, "plots", "ensembl_orthologs.pdf"),
    interaction_plot=join(PROCESSED_DIR, "plots", "cellphonedb_interactions.pdf")
  notebook:
    join("src", "cellphonedb_orthologs.py.ipynb")

rule metmap_pca:
  input:
    mm_potential=join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx")
  output:
    pca_plot=join(PROCESSED_DIR, "metmap", "pca.pdf")
  notebook:
    join("src", "metmap_pca.py.ipynb")

rule convert_gct_to_h5ad:
  input:
    join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.gct")
  output:
    join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.h5ad")
  script:
    join("src", "convert_gct_to_h5ad.py")
    
rule download_orthologs:
  output:
    join(RAW_DIR, "ensembl", "human_mouse_orthologs.tsv")
  script:
    join("src", "download_orthologs.py")

rule generate_pseudobulk:
  input:
    join(RAW_DIR, "tm", "seurat", "{tissue}.facs.Robj")
  output:
    join(RAW_DIR, "tm", "anndata", "{tissue}.pseudobulk.h5ad")
  script:
    join("src", "generate_pseudobulk.R")

rule convert_robj_to_h5ad:
  input:
    join(RAW_DIR, "tm", "seurat", "{tissue}.facs.Robj")
  output:
    join(RAW_DIR, "tm", "anndata", "{tissue}.facs.h5ad")
  script:
    join("src", "convert_robj_to_h5ad.R")


# Abstract parent rules
rule curl_download:
  shell:
    '''
    curl -L -o {output} "{params.file_url}"
    '''

rule unzip:
  shell:
    '''
    unzip -o {input} -d {params.out_dir}
    '''

rule gunzip:
  shell:
    '''
    gunzip -c {input} > {output}
    '''
    
# Download Tabula Muris data
use rule curl_download as download_tm_facs_robj with:
  output:
    join(RAW_DIR, "tm", "seurat", "{tissue}.facs.Robj")
  params:
    file_url=(lambda w: TM_FACS_SEURAT[w.tissue])

# Download CellPhoneDB data
use rule curl_download as download_cellphonedb_gene_input with:
  output:
    join(RAW_DIR, "cellphonedb", "gene_input.csv")
  params:
    file_url=CELLPHONEDB_GENE_INPUT_URL

use rule curl_download as download_cellphonedb_protein_input with:
  output:
    join(RAW_DIR, "cellphonedb", "protein_input.csv")
  params:
    file_url=CELLPHONEDB_PROTEIN_INPUT_URL
  
use rule curl_download as download_cellphonedb_complex_input with:
  output:
    join(RAW_DIR, "cellphonedb", "complex_input.csv")
  params:
    file_url=CELLPHONEDB_COMPLEX_INPUT_URL

use rule curl_download as download_cellphonedb_interaction_input with:
  output:
    join(RAW_DIR, "cellphonedb", "interaction_input.csv")
  params:
    file_url=CELLPHONEDB_INTERACTION_INPUT_URL
    
# Download MetMap data
use rule curl_download as download_metmap_500 with:
  output:
    join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx")
  params:
    file_url=METMAP_500_URL

use rule gunzip as gunzip_metmap_counts with:
  input:
    join(RAW_DIR, "metmap", "GSE148283_all.count.csv.gz")
  output:
    join(RAW_DIR, "metmap", "GSE148283_all.count.csv")

use rule gunzip as gunzip_metmap_samples with:
  input:
    join(RAW_DIR, "metmap", "GSE148283_all.sample.csv.gz")
  output:
    join(RAW_DIR, "metmap", "GSE148283_all.sample.csv")

use rule curl_download as download_metmap_counts with:
  output:
    join(RAW_DIR, "metmap", "GSE148283_all.count.csv.gz")
  params:
    file_url=METMAP_EXP_COUNTS_URL

use rule curl_download as download_metmap_samples with:
  output:
    join(RAW_DIR, "metmap", "GSE148283_all.sample.csv.gz")
  params:
    file_url=METMAP_EXP_SAMPLES_URL
    
# Download CCLE data
use rule gunzip as gunzip_ccle_exp with:
  input:
    join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.gct.gz")
  output:
    join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.gct")
    
use rule curl_download as download_ccle_exp with:
  output:
    join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.gct.gz")
  params:
    file_url=CCLE_EXP_URL

# Download ontology files
use rule curl_download as download_cl_obo with:
  output:
    join(RAW_DIR, "ontologies", "cl.obo")
  params:
    file_url=CL_OBO_URL
  
use rule curl_download as download_go_obo with:
  output:
    join(RAW_DIR, "ontologies", "go.obo")
  params:
    file_url=GO_OBO_URL