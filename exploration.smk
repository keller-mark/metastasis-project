
include: "downloads.smk"

configfile: "config.yml"

assert("params" in config.keys())

PARAMS = config["params"]

# Abbreviations
# - iea: interaction_expression_aggregation (how to combine expression for partners A and B of an interaction)
# - cea: complex_expression_aggregation (how to aggregate expression for complexes which contain multiple gene products)
# - fic: feature_inclusion_criteria (which features to include when building the model)

# Rules
rule all:
  input:
    expand(
      join(PROCESSED_DIR, "nonmet", "{tissue}.gsea.pos.tsv"),
      tissue=METMAP_TISSUES,
    ),
    
# Run differential expression and gene set enrichment methods
# on non-metastatic vs. metastatic samples for each target organ type.
rule nonmetastatic_comparison:
  input:
    deseq=join(PROCESSED_DIR, "nonmet", "{tissue}.deseq.results.csv"),
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue])
  output:
    deseq_plot=join(PROCESSED_DIR, "nonmet", "{tissue}.deseq.plot.pdf"),
    enrichr=join(PROCESSED_DIR, "nonmet", "{tissue}.enrichr.tsv"),
    enrichr_plot=join(PROCESSED_DIR, "nonmet", "{tissue}.enrichr.plot.pdf"),
    enrichr_pos_plot=join(PROCESSED_DIR, "nonmet", "{tissue}.enrichr.pos_plot.pdf"),
    enrichr_neg_plot=join(PROCESSED_DIR, "nonmet", "{tissue}.enrichr.neg_plot.pdf"),
    gsea_pos=join(PROCESSED_DIR, "nonmet", "{tissue}.gsea.pos.tsv"),
    gsea_neg=join(PROCESSED_DIR, "nonmet", "{tissue}.gsea.neg.tsv"),
    gsea_pos_plot=join(PROCESSED_DIR, "nonmet", "{tissue}.gsea.pos_plot.pdf"),
    gsea_neg_plot=join(PROCESSED_DIR, "nonmet", "{tissue}.gsea.neg_plot.pdf"),
  notebook:
    join("src", "nonmetastatic_comparison.py.ipynb")

rule nonmetastatic_deseq:
  input:
    counts=join(PROCESSED_DIR, "nonmet", "{tissue}.deseq.counts.csv"),
    conditions=join(PROCESSED_DIR, "nonmet", "{tissue}.deseq.conditions.csv"),
  output:
    join(PROCESSED_DIR, "nonmet", "{tissue}.deseq.results.csv"),
  script:
    join("src", "nonmetastatic_deseq.R")

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
rule nonmetastatic_deseq_preparation:
  input:
    mm_potential=join(RAW_DIR, "metmap", "metmap_500_met_potential.xlsx"),
    ccle_exp=join(RAW_DIR, "ccle", "CCLE_RNAseq_genes_counts_20180929.h5ad"),
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue])
  output:
    counts=join(PROCESSED_DIR, "nonmet", "{tissue}.deseq.counts.csv"),
    conditions=join(PROCESSED_DIR, "nonmet", "{tissue}.deseq.conditions.csv"),
  notebook:
    join("src", "nonmetastatic_deseq_preparation.py.ipynb")

# Use co-expression values to build a model of metastasis potential
# for each MetMap target site.
rule build_model:
  input:
    join(PROCESSED_DIR, "coexpression", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression.h5ad"),
  params:
    metmap_tissue=(lambda w: TM_TO_METMAP[w.tissue])
  output:
    model=join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{fic}.{tissue}.model.h5ad"),
    mse_plot=join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{fic}.{tissue}.mse_plot.pdf"),
    prediction_plot=join(PROCESSED_DIR, "models", "{model}.{expression_scale}.{interaction_source}.{cea}.{iea}.{fic}.{tissue}.prediction_plot.pdf")
  notebook:
    join("src", "build_model.py.ipynb")
  
# Plot distribution of Gini coefficients across interaction coexpression values for each cancer cell line, for each tissue.
rule plot_coexpression_variation:
  input:
    coexpression=join(PROCESSED_DIR, "coexpression", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression.h5ad"),
    interactions=join(PROCESSED_DIR, "cellphonedb", "gene_orthologs.tsv")
  params:
    tissues=METMAP_TISSUES
  output:
    gini=join(PROCESSED_DIR, "plots", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression_gini.pdf"),
    entropy=join(PROCESSED_DIR, "plots", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression_entropy.pdf"),
    variance=join(PROCESSED_DIR, "plots", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression_variance.pdf"),
  notebook:
    join("src", "plot_coexpression_variation.py.ipynb")
    
# Plot top co-expression values for each (cancer cell line, tabula muris cell type) pair.
rule plot_top_coexpression:
  input:
    coexpression=join(PROCESSED_DIR, "coexpression", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression.h5ad"),
    interactions=join(PROCESSED_DIR, "cellphonedb", "gene_orthologs.tsv")
  params:
    tissues=METMAP_TISSUES
  output:
    join(PROCESSED_DIR, "plots", "{expression_scale}.{interaction_source}.{cea}.{iea}.{tissue}.coexpression_top_5.html")
  notebook:
    join("src", "plot_top_coexpression.py.ipynb")

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
    pca_plot_1=join(PROCESSED_DIR, "plots", "metmap_pca_1.pdf"),
    pca_plot_2=join(PROCESSED_DIR, "plots", "metmap_pca_2.pdf"),
    pca_plot_3=join(PROCESSED_DIR, "plots", "metmap_pca_3.pdf"),
    pca_plot_4=join(PROCESSED_DIR, "plots", "metmap_pca_4.pdf"),
    pca_plot_5=join(PROCESSED_DIR, "plots", "metmap_pca_5.pdf"),
    pca_plot_6=join(PROCESSED_DIR, "plots", "metmap_pca_6.pdf")
  notebook:
    join("src", "metmap_pca.py.ipynb")

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
