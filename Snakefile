from os.path import join
from urllib.parse import quote_plus

DATA_DIR = "data"
RAW_DIR = join(DATA_DIR, "raw")
PROCESSED_DIR = join(DATA_DIR, "processed")

TM_MATRICES_URL = "https://storage.googleapis.com/broad-datarepo-terra-prod-hca2-bucket/d30e68f8-c826-4639-88f3-ae35f00d4185/2f393f0c-3ff2-4ce3-9fa4-b4bcf97b3a96/e0009214-c0a0-4a7b-96e2-d6a83e966ce0.FACS_matrices.zip?Expires=1614305822&GoogleAccessId=azul-ucsc-0-prod%40platform-hca-prod.iam.gserviceaccount.com&Signature=XP8uCE81Ly6hI4TfUbsSdANyh%2BDkMlZEa1OW9mGV8pV0q6X2icjby6VnrMsJdjVCrV1HQyOmJOKZ7ID1QY%2Bv2v9Hq7ZIkbNX92FSkoCKHyFjh0kIPQvpmXbbF7feuPpWqtfipPKjk6u3JvTe2wzJu11TxqrkjXXl0fYlpQkhV%2FJ1cX%2FCELGetChSY4m7QfrW8Mhwoh8YVxwmiZxbZP8Bo9jgz4vS6iGKRb%2F98O1Aa3eBgSFLL26zu2JFeMFtTrzSM3y04u%2FQD7JHSqL4m89LqmgjxkE96KkA5c0X3PtrXFLKdsXx%2BAl9xLOvacgqmL3mRBYyMdJdqD2H8ipFoVUi9w%3D%3D&response-content-disposition=attachment%3B+filename%3De0009214-c0a0-4a7b-96e2-d6a83e966ce0.FACS_matrices.zip"

TM_ANNOTATIONS_URL = "https://storage.googleapis.com/broad-datarepo-terra-prod-hca2-bucket/d30e68f8-c826-4639-88f3-ae35f00d4185/2876118f-3607-4718-8f82-7a53fa0ecef4/e0009214-c0a0-4a7b-96e2-d6a83e966ce0.annotations_facs.csv?Expires=1614295616&GoogleAccessId=azul-ucsc-0-prod%40platform-hca-prod.iam.gserviceaccount.com&Signature=VDlgr6gnITqGI5nGR0FW5Kt4jirZ61w4Th5CR%2B19eCUXgoJQYRmpMnakXJCqXB1iZcawnY1bm7Scw5MTrXiX25pIN8vH5xGrEjNcyDaeFICHGlZ%2B3nHdg0jjHLy9PNaIs9tLeY7uyMGmpdfBPpc%2BmDX7JPce4AJjIC6xD4DOULdgGHaBDz4YqQLtt7V5pw13h2bXFQZSfeX3HPH2oKXGqOa80eaAbZAyDhYp8S8H8bVPRmtQnnpGCfLVy8zZ%2F2mQjZr8kOlnsDSJQVSr0deawFEuXxZ93mk76%2FW6Tcg3LC1APXpeC%2FK4E9tu4eyQV7Q5SICYDG%2Fs6HPrJ5T0XcsugA%3D%3D&response-content-disposition=attachment%3B+filename%3De0009214-c0a0-4a7b-96e2-d6a83e966ce0.annotations_facs.csv"

  
rule all:
  input:
    join(RAW_DIR, "tabula_muris", "FACS", "Aorta-counts.csv"),
    join(RAW_DIR, "tabula_muris", "annotations.csv")

rule download_annotations:
  output:
      join(RAW_DIR, "tabula_muris", "annotations.csv")
  params:
      file_url=TM_MATRICES_URL
  shell:
    '''
    curl -L -o {output} "{params.file_url}"
    '''

rule unzip_matrices:
  input:
    join(RAW_DIR, "tabula_muris", "matrices.zip")
  output:
    join(RAW_DIR, "tabula_muris", "FACS", "Aorta-counts.csv")
  shell:
    '''
    unzip {input} -d $(dirname $(dirname {output}))
    '''
    


rule download_matrices:
  output:
      join(RAW_DIR, "tabula_muris", "matrices.zip")
  params:
      file_url=TM_MATRICES_URL
  shell:
    '''
    curl -L -o {output} "{params.file_url}"
    '''