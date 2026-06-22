# Paired-diagnosis-relapse-pAML-data-analysis
This repository contains the custom R code used to generate the analyses presented in the manuscript:  "Characterization of Chemoresistant Cell Populations Improves Risk Stratification and Therapy Prediction in Pediatric Acute Myeloid Leukemia"
The repository was simplified from the original analysis environment to provide a reproducible version of the major analyses used in the study.

The code is organized into three primary modules:

Single-cell RNA-seq preprocessing and integration
pAML subpopulation identification
Downstream validation and survival analyses

Previously published methods (e.g. SQUID deconvolution) are not reimplemented here and should be obtained from their original publication/repository.

Data requirements

This repository assumes access to the datasets used in the manuscript.

Input data may include:

Single-cell count matrices
Patient metadata
External validation cohorts (e.g. TARGET AML)
Clinical outcome tables

As the datasets are large, they are not included in this repository.

Required R packages
Seurat
SingleR

dplyr
tidyr
data.table

ggplot2
stringr

survival
survminer

Install packages before running analyses.

Citation

If using this code, please cite: NajafPanah MJ, Stevens AM, Krueger MJ, Rochette M, Sandhu S, Kim L, Chiu HS, Epps J, Somvanshi S, Zorman B, Martinez MR, Rapsomaniki M, Unger S, Becher B, Yi JS, Man TK, Redell MS, Sumazin P. Characterization of Chemoresistant Cell Populations Improves Risk Stratification and Therapy Prediction in Pediatric AML. bioRxiv [Preprint]. 2026 Oct 20:2025.09.25.678688. doi: 10.1101/2025.09.25.678688. PMID: 41256415; PMCID: PMC12622046. In Press, Nature Communications.
