# MotrpacRatTraining6moQCRep

### About this repository 
This repository provides R Markdown scripts (.Rmd) and corresponding HTML reports with 
the post-quantification analyses of the data presented in the main paper for the first 
large-scale multi-omic multi-tissue endurance exercise training study conducted 
in young adult rats by the Molecular Transducers of Physical Activity Consortium 
(MoTrPAC). Find the [preprint on bioRxiv](https://doi.org/10.1101/2022.09.21.508770).

Each omic data generated for the paper above had a pipeline to generate quantified features from the raw data (e.g., gene-level counts for transcriptomics data quantified from raw RNA-seq reads). In this repository we provide scripts with the initial analyses of these post-quantification data, which are exploratory in nature. These analyses were used for sanity checks, QA procedures, and identifying flagged samples such as extreme PCA outliers. 

The reports provided in this repository were originally intended for use within MoTrPAC. Therefore, some functionalities that fetch data directly from the internal MoTrPAC cloud systems are not expected to work as-is. Nevertheless, we decided to provide these reports because: (1) they contain useful code/analyses and figures (e.g., PC plots of each tissue in the [transcriptomics report](https://htmlpreview.github.io/?https://github.com/MoTrPAC/MotrpacRatTraining6moQCRep/blob/main/reports/pass1b_rnaseq.html)), and (2) most analyses can be reproduced by downloading the data from the [MoTrPAC Data Hub](https://motrpac-data.org/) and modifying the code to take the data from the local directories to which the data were downloaded. The relevant code sections and paths in each report (and its Rmd file) typically use `ADDPATH` as a token that replaced internal MoTrPAC cloud paths.

### About MoTrPAC
MoTrPAC is a national research consortium designed to discover and perform 
preliminary characterization of the range of molecular transducers (the 
"molecular map") that underlie the health benefits of physical activity in humans. 
The program's goal is to study the molecular changes that occur during and after 
exercise and ultimately to advance the understanding of how physical activity 
improves and preserves health. The six-year program is the largest targeted NIH 
investment of funds into the mechanisms of how physical activity improves health 
and prevents disease. See [motrpac.org](https://www.motrpac.org/) and 
[motrpac-data.org](https://motrpac-data.org/) for more details. 

### Repository structure

This repo is organized as follows:

1. `tools` - functions, can source each `.R` file into a session
1. `reports` - the main folder with the reports (`.Rmd` files and the corresponding HTML reports)
1. `ext_data` - additional data required for reproducing some analyses (e.g., genome.gtf for the ATAC-seq analysis)

### Links to rendered HTML reports 

**NOTE:** These links will work when the repository is made public. 

The HTML reports can be viewed either by downloading this repository and opening the HTML files or by clicking the following links.

- [Chromatin accessibilty (ATAC-seq)](https://htmlpreview.github.io/?https://github.com/MoTrPAC/MotrpacRatTraining6moQCRep/blob/main/reports/pass1b_atacseq.html)
- [Multiplexed immunoassays](https://htmlpreview.github.io/?https://github.com/MoTrPAC/MotrpacRatTraining6moQCRep/blob/main/reports/pass1b_immunoassay.html)
- [Proteomics](https://htmlpreview.github.io/?https://github.com/MoTrPAC/MotrpacRatTraining6moQCRep/blob/main/reports/pass1b_proteomics.html)
- [Transcriptomics (RNA-seq)](https://htmlpreview.github.io/?https://github.com/MoTrPAC/MotrpacRatTraining6moQCRep/blob/main/reports/pass1b_rnaseq.html)
- [DNA Methylation (RRBS)](https://htmlpreview.github.io/?https://github.com/MoTrPAC/MotrpacRatTraining6moQCRep/blob/main/reports/pass1b_rrbs.html)
