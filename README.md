Prediction of metastasis target site using ligand-receptor co-expression (human cancer cell lines in vitro and healthy organs in vivo)

- How accurately can metastasis target sites be predicted given co-expression of ligands and receptors between human cancer cell lines and healthy mouse tissues?
  - Can a "lower bound" be placed on the seed-and-soil hypothesis (Paget 1889)?
    - Obviously, there are many other factors to consider besides the ligand-receptor interactions between circulating cancer cells and healthy tissues. If a more complex model is able to take into account factors such as physical forces, vasculature (blood and lymphatic), metabolic differences, etc. how much does each factor improve prediction accuracy?

[PDF](https://github.com/keller-mark/lr/blob/gh-pages/main.pdf), [Notes](./NOTES.md)

## Set up environment

```sh
conda env create -f environment.yml
conda activate lr-env
```

```R
install.packages("BiocManager")
install.packages("devtools")
install.packages("reticulate")
install.packages("Seurat")
install.packages("jsonlite")
BiocManager::install("splatter")
BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))
devtools::install_github("HelenaLC/muscat", ref = "master")
devtools::install_github("cellgeni/sceasy")
```

## Run pipeline

```sh
snakemake -j 1
```

## Edit a notebook in the pipeline

```sh
snakemake -j 1 --edit-notebook data/intermediate/coexpression/Kidney.coexpression.h5ad
```

## Resources

- [Tabula muris](https://tabula-muris.ds.czbiohub.org/): 20 organs and tissues from mouse
  - [Human Cell Atlas data portal](https://data.humancellatlas.org/explore/projects/e0009214-c0a0-4a7b-96e2-d6a83e966ce0/expression-matrices?catalog=dcp2)
  - Zhang et al. Cell 2013: "The evidence supports a model in which a TN tumor stroma rich in mesenchymal cytokines CXCL12 and IGF1 skews the carcinoma population toward a preponderance of clones with a high level of activated Src and a predisposition to grow in the bone marrow, a soil that is rich in CXCL12 and IGF1 compared to lung, brain, and other metastasis target tissues"
    - task: explore tabula muris to determine whether mouse bone marrow data seems to have this property as well

- [MetMap downloads](https://depmap.org/metmap/data/index.html)
- [Cancer cell line encyclopedia (CCLE)](https://portals.broadinstitute.org/ccle)
- [CellPhoneDB](https://www.cellphonedb.org/)
