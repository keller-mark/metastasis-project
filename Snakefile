from os.path import join
from urllib.parse import quote_plus

DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
PROCESSED_DIR = join(DATA_DIR, "processed")

# Tabula Muris URLs: https://data.humancellatlas.org/explore/projects/e0009214-c0a0-4a7b-96e2-d6a83e966ce0/expression-matrices?catalog=dcp2
TM_MATRICES_URL = "https://storage.googleapis.com/broad-datarepo-terra-prod-hca2-bucket/d30e68f8-c826-4639-88f3-ae35f00d4185/2f393f0c-3ff2-4ce3-9fa4-b4bcf97b3a96/e0009214-c0a0-4a7b-96e2-d6a83e966ce0.FACS_matrices.zip?Expires=1614305822&GoogleAccessId=azul-ucsc-0-prod%40platform-hca-prod.iam.gserviceaccount.com&Signature=XP8uCE81Ly6hI4TfUbsSdANyh%2BDkMlZEa1OW9mGV8pV0q6X2icjby6VnrMsJdjVCrV1HQyOmJOKZ7ID1QY%2Bv2v9Hq7ZIkbNX92FSkoCKHyFjh0kIPQvpmXbbF7feuPpWqtfipPKjk6u3JvTe2wzJu11TxqrkjXXl0fYlpQkhV%2FJ1cX%2FCELGetChSY4m7QfrW8Mhwoh8YVxwmiZxbZP8Bo9jgz4vS6iGKRb%2F98O1Aa3eBgSFLL26zu2JFeMFtTrzSM3y04u%2FQD7JHSqL4m89LqmgjxkE96KkA5c0X3PtrXFLKdsXx%2BAl9xLOvacgqmL3mRBYyMdJdqD2H8ipFoVUi9w%3D%3D&response-content-disposition=attachment%3B+filename%3De0009214-c0a0-4a7b-96e2-d6a83e966ce0.FACS_matrices.zip"

TM_ANNOTATIONS_URL = "https://storage.googleapis.com/broad-datarepo-terra-prod-hca2-bucket/d30e68f8-c826-4639-88f3-ae35f00d4185/2876118f-3607-4718-8f82-7a53fa0ecef4/e0009214-c0a0-4a7b-96e2-d6a83e966ce0.annotations_facs.csv?Expires=1614307121&GoogleAccessId=azul-ucsc-0-prod%40platform-hca-prod.iam.gserviceaccount.com&Signature=fSAhMP06zmH2zmIxp6QI70%2FVuYon2by8uip6Al956Z7j2tquc8QboOqKfL2lU74aHLXcdHoZPmepvQ1ozrNdGwaozokFotkQ2OYUsFna3VwUHz11jFgPlk7n5kUyWhykRRS8NgtfGr21jmH3QghwKNBAban%2F%2FeZTRK32N18qr4ap0%2FvCYfGRtBGO4WH%2FbgvjPVVxvhOr8Jq4fxYfNXK8iaiUy71pXoFjPJilq3PXOAw7cAlkPNuav6n9I5D2NIXUy8brEoZoMI%2BJHicVdVwJMOP%2FO2rwPQPHs6vJGLRFIjfwwAKWu9vwOi%2FG%2BfMCq4PtDFeONDpgCvcq%2BCzt5exIrA%3D%3D&response-content-disposition=attachment%3B+filename%3De0009214-c0a0-4a7b-96e2-d6a83e966ce0.annotations_facs.csv"

TM_TISSUES = [
 'Bladder',
 'Brain_Non-Myeloid',
 'Fat',
 'Skin',
 'Trachea',
 'Lung',
 'Kidney',
 'Brain_Myeloid',
 'Pancreas',
 'Tongue',
 'Limb_Muscle',
 'Thymus',
 'Mammary_Gland',
 'Diaphragm',
 'Marrow',
 'Spleen',
 'Large_Intestine',
 'Liver',
 'Heart',
 'Aorta'
]

# CellPhoneDB URLs: https://www.cellphonedb.org/downloads
CELLPHONEDB_GENE_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv"
CELLPHONEDB_PROTEIN_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/protein_input.csv"
CELLPHONEDB_PROTEIN_CURATED_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/protein_curated.csv"

# MetMap URLs: https://depmap.org/metmap/data/index.html
METMAP_500_URL = "https://ndownloader.figshare.com/files/24009293"

# Rules
rule all:
  input:
    expand(join(RAW_DIR, "tm", "adata", "{tissue}.h5ad"), tissue=TM_TISSUES),
    join(RAW_DIR, "cellphonedb", "gene_input.csv"),
    join(RAW_DIR, "cellphonedb", "protein_input.csv"),
    join(RAW_DIR, "cellphonedb", "protein_curated.csv"),
    join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx")
    

rule convert_to_h5ad:
  input:
    counts=join(RAW_DIR, "tm", "FACS", "{tissue}-counts.csv"),
    annots=join(RAW_DIR, "tm", "annotations.csv")
  output:
    join(RAW_DIR, "tm", "adata", "{tissue}.h5ad")
  script:
    join("src", "convert_to_h5ad.py")



rule unzip_matrices:
  input:
    join(RAW_DIR, "tm", "matrices.zip")
  output:
    expand(join(RAW_DIR, "tm", "FACS", "{tissue}-counts.csv"), tissue=TM_TISSUES)
  params:
    out_dir=join(RAW_DIR, "tm")
  shell:
    '''
    unzip -o {input} -d {params.out_dir}
    '''

# Abstract parent rules
rule curl_download:
  shell:
    '''
    curl -L -o {output} "{params.file_url}"
    '''

# Download Tabula Muris data
use rule curl_download as download_tm_matrices with:
  output:
    join(RAW_DIR, "tm", "matrices.zip")
  params:
    file_url=TM_MATRICES_URL

use rule curl_download as download_tm_annotations with:
  output:
    join(RAW_DIR, "tm", "annotations.csv")
  params:
    file_url=TM_ANNOTATIONS_URL

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
  
use rule curl_download as download_cellphonedb_protein_curated with:
  output:
    join(RAW_DIR, "cellphonedb", "protein_curated.csv")
  params:
    file_url=CELLPHONEDB_PROTEIN_CURATED_URL
    
# Download MetMap data
use rule curl_download as download_metmap_500 with:
  output:
    join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx")
  params:
    file_url=METMAP_500_URL