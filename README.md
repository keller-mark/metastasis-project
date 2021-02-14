# metastasis-prediction-by-ligand-receptor

Goal: determine whether ligand-receptor interactions are predictive for metastasis destination organ or tissue type.

Ground truth:
- for each cancer type, observed frequencies of each organ or tissue type as the metastasis destination

Input data:
- scRNA-seq data for each primary cancer type
- scRNA-seq data for each healthy organ or tissue type
- ligand-receptor pairs

Pre-processing:
- identify the immune/stromal cell types environment in each healthy organ
- identify the ligand and receptor signals in the healthy immune cells
- identify the ligand and receptor signals in the cancer samples


Model inputs: for each ligand and receptor pair, abundance of RNA in cancer cells vs each healthy organ

Model output: score for metastasis be localized to each healthy organ (highest score implies most likely for metastasis destination)

Questions:
- Could a graph embedding be used?
  - Perhaps for each cancer type, build an input graph of ligand-receptor pairs (edges between pairs across organs, with weights for interaction strength), embed to map to a graph of cancer type to organ type with weights
