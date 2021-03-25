
include: "downloads.smk"

NUM_FOLDS = 5

# Want to do 5-fold cross validation to obtain differentially expressed genes and pathways
# in each fold, then train model with the data in the fold and predict on the held-out fold

rule all:
  input:
    expand(
      join(PROCESSED_DIR, "kfold_deseq", "{fold}.deseq.model.json"),
      fold=range(NUM_FOLDS),
    )

# Build a PLSRegression model for each fold, using the union of significantly
# differentially expressed genes across tissue types (discovered in the training set).
rule kfold_deseq_model:
  input:
    deseq=expand(
      join(PROCESSED_DIR, "kfold_deseq", "{{fold}}.{tissue}.deseq.results.csv"),
      tissue=METMAP_TISSUES
    ),
    kfold_indices=join(PROCESSED_DIR, "kfold_deseq", "kfold_indices.csv"),
    mm_potential=join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx"),
    ccle_exp=join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.h5ad")
  params:
    metmap_tissues=METMAP_TISSUES,
    tm_to_metmap=TM_TO_METMAP,
  output:
    model_plot=join(PROCESSED_DIR, "kfold_deseq", "{fold}.deseq.model.pdf"),
    model_results=join(PROCESSED_DIR, "kfold_deseq", "{fold}.deseq.model.json"),
  notebook:
    join("src", "kfold_deseq_model.py.ipynb")


# Plot differential expression results for each fold and tissue type.
rule kfold_deseq_plot:
  input:
    deseq=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.results.csv")
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue]),
  output:
    deseq_plot=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.plot.pdf"),
    enrichr=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.enrichr.tsv"),
    enrichr_plot=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.enrichr.plot.pdf")
  notebook:
    join("src", "kfold_deseq_plot.py.ipynb")

# Run differential expression and gene set enrichment methods
# on non-metastatic vs. metastatic samples for each target organ type.
rule kfold_deseq:
  input:
    counts=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.counts.csv"),
    conditions=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.conditions.csv"),
  output:
    join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.results.csv"),
  script:
    join("src", "kfold_deseq.R")

# Run differential expression on the 4/5 of the data that is the training set.
# Hold out the remaining 1 fold for testing.
rule kfold_deseq_preparation:
  input:
    kfold_indices=join(PROCESSED_DIR, "kfold_deseq", "kfold_indices.csv"),
    mm_potential=join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx"),
    ccle_exp=join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.h5ad"),
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue]),
  output:
    counts=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.counts.csv"),
    conditions=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.conditions.csv"),
  notebook:
    join("src", "kfold_deseq_preparation.py.ipynb")

# Run differential expression on the 4/5 of the data that is the training set.
# Hold out the remaining 1 fold for testing.
rule kfold_indices:
  input:
    mm_potential=join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx"),
    ccle_exp=join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.h5ad"),
  params:
    num_folds=NUM_FOLDS,
  output:
    join(PROCESSED_DIR, "kfold_deseq", "kfold_indices.csv"),
  notebook:
    join("src", "kfold_indices.py.ipynb")