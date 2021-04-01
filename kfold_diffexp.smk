
include: "downloads.smk"

NUM_FOLDS = 5

DESEQ_FC_THRESHOLD = [
  2
]

NUM_PCS = [
  5
]

GEXP_TRANSFORMS = [
  "tpm",
  "log1p_tpm"
]

# Want to do 5-fold cross validation to obtain differentially expressed genes and pathways
# in each fold, then train model with the data in the fold and predict on the held-out fold

rule all:
  input:
    expand(
      join(PROCESSED_DIR, "kfold_deseq", "{num_pc}.{gexp_transform}.{fc_threshold}.deseq.model.performance.pdf"),
      num_pc=NUM_PCS,
      gexp_transform=GEXP_TRANSFORMS,
      fc_threshold=DESEQ_FC_THRESHOLD
    ),
    expand(
      join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.heatmap.plot.pdf"),
      fold=range(NUM_FOLDS),
      tissue=METMAP_TISSUES,
    )

rule kfold_deseq_model_plot:
  input:
    expand(
      join(PROCESSED_DIR, "kfold_deseq", "{fold}.{{num_pc}}.{{gexp_transform}}.{{fc_threshold}}.deseq.model.json"),
      fold=range(NUM_FOLDS),
    )
  params:
    metmap_tissues=METMAP_TISSUES,
    tm_to_metmap=TM_TO_METMAP,
  output:
    plot=join(PROCESSED_DIR, "kfold_deseq", "{num_pc}.{gexp_transform}.{fc_threshold}.deseq.model.performance.pdf"),
  notebook:
    join("src", "kfold_deseq_model_plot.py.ipynb")

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
    model_results=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{num_pc}.{gexp_transform}.{fc_threshold}.deseq.model.json"),
  notebook:
    join("src", "kfold_deseq_model.py.ipynb")

rule kfold_deseq_heatmap:
  input:
    counts=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.counts.csv"),
    conditions=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.conditions.csv"),
    deseq=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.deseq.results.csv")
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue]),
  output:
    heatmap_plot=join(PROCESSED_DIR, "kfold_deseq", "{fold}.{tissue}.heatmap.plot.pdf")
  notebook:
    join("src", "kfold_deseq_heatmap.py.ipynb")
    
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