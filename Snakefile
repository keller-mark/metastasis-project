from os.path import join
from urllib.parse import quote_plus

DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
INTERMEDIATE_DIR = join(DATA_DIR, "intermediate")
PROCESSED_DIR = join(DATA_DIR, "processed")

TM_FACS_ZIP_URL = "https://ndownloader.figshare.com/files/10700143"
TM_FACS_ANNOTATIONS_URL = "https://ndownloader.figshare.com/files/13088129"
TM_FACS_METADATA_URL = "https://ndownloader.figshare.com/files/10842785"

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

# CellPhoneDB URLs: https://www.cellphonedb.org/downloads
CELLPHONEDB_GENE_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/gene_input.csv"
CELLPHONEDB_PROTEIN_INPUT_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/protein_input.csv"
CELLPHONEDB_PROTEIN_CURATED_URL = "https://raw.githubusercontent.com/Teichlab/cellphonedb-data/master/data/sources/protein_curated.csv"

# MetMap URLs: https://depmap.org/metmap/data/index.html
METMAP_500_URL = "https://ndownloader.figshare.com/files/24009293"

# Rules
rule all:
  input:
    expand(join(RAW_DIR, "tm", "seurat", "{tissue}.facs.Robj"), tissue=TM_FACS_TISSUES),
    expand(join(INTERMEDIATE_DIR, "tm", "splatter", "{tissue}.params.json"), tissue=TM_FACS_TISSUES),
    expand(join(RAW_DIR, "tm", "anndata", "{tissue}.facs.h5ad"), tissue=TM_FACS_TISSUES),
    join(RAW_DIR, "cellphonedb", "gene_input.csv"),
    join(RAW_DIR, "cellphonedb", "protein_input.csv"),
    join(RAW_DIR, "cellphonedb", "protein_curated.csv"),
    join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx")

rule estimate_splatter_params:
  input:
    join(RAW_DIR, "tm", "seurat", "{tissue}.facs.Robj")
  output:
    join(INTERMEDIATE_DIR, "tm", "splatter", "{tissue}.params.json")
  script:
    join("src", "estimate_splatter_params.R")

rule sample_bulk_from_tm_sc_by_tissue:
  input:
    join(RAW_DIR, "tm", "anndata", "{tissue}.facs.h5ad")
  output:
    join(INTERMEDIATE_DIR, "tm", "anndata", "{tissue}.bulk.h5ad")
  script:
    join("src", "sample_bulk_from_tm_sc_by_tissue.py")

rule convert_tm_to_h5ad:
  input:
    counts=join(RAW_DIR, "tm", "FACS", "{tissue}-counts.csv"),
    annots=join(RAW_DIR, "tm", "facs_annotations.csv")
  output:
    join(RAW_DIR, "tm", "anndata", "{tissue}.facs.h5ad")
  script:
    join("src", "convert_tm_to_h5ad.py")


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

# Download Tabula Muris data
use rule curl_download as download_tm_facs_robj with:
  output:
    join(RAW_DIR, "tm", "seurat", "{tissue}.facs.Robj")
  params:
    file_url=(lambda w: TM_FACS_SEURAT[w.tissue])

use rule curl_download as download_tm_facs_zip with:
  output:
    join(RAW_DIR, "tm", "facs_matrices.zip")
  params:
    file_url=TM_FACS_ZIP_URL

use rule unzip as unzip_tm_facs_zip with:
  input:
    join(RAW_DIR, "tm", "facs_matrices.zip")
  output:
    expand(join(RAW_DIR, "tm", "FACS", "{tissue}-counts.csv"), tissue=TM_FACS_TISSUES)
  params:
    out_dir=join(RAW_DIR, "tm")

use rule curl_download as download_tm_facs_annotations with:
  output:
    join(RAW_DIR, "tm", "facs_annotations.csv")
  params:
    file_url=TM_FACS_ANNOTATIONS_URL

use rule curl_download as download_tm_facs_metadata with:
  output:
    join(RAW_DIR, "tm", "facs_metadata.csv")
  params:
    file_url=TM_FACS_METADATA_URL

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