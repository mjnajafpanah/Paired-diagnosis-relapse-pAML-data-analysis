# Paired-diagnosis-relapse-pAML-data-analysis
This repository contains the custom R code used to generate the analyses presented in "Characterization of Chemoresistant Cell Populations Improves Risk Stratification and Therapy Prediction in Pediatric Acute Myeloid Leukemia" research paper.
<p align="center">
  <img width="424" height="420" alt="Image" src="https://github.com/user-attachments/assets/0c5bc333-fa22-469c-8aa8-e92b89ff275f" />
</p>

The repository was simplified from the original analysis environment to provide a reproducible version of the major analyses used in the study. The code is organized into three primary modules:

1. Single-cell RNA-seq preprocessing and integration
2. pAML subpopulation identification
3. Downstream validation and survival analyses

Previously published methods (e.g. SQUID deconvolution) are not reimplemented here and should be obtained from their original publication/repository.

## Input data

- Single-cell count matrices
- Patient metadata
- External validation cohorts (e.g. TARGET AML)
- Clinical outcome tables

As the datasets are large, they are not included in this repository.

## Required R packages
- Seurat
- SingleR
- dplyr
- tidyr
- data.table
- ggplot2
- stringr
- survival
- survminer

Install packages before running analyses.

## Citation

If using this code, please cite: 
NajafPanah MJ, Stevens AM, Krueger MJ, Rochette M, Sandhu S, Kim L, Chiu HS, Epps J, Somvanshi S, Zorman B, Martinez MR, Rapsomaniki M, Unger S, Becher B, Yi JS, Man TK, Redell MS, Sumazin P. Characterization of Chemoresistant Cell Populations Improves Risk Stratification and Therapy Prediction in Pediatric AML. bioRxiv (2025). doi:10.1101/2025.09.25.678688. Accepted for publication in Nature Communications.
