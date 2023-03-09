# MotrpacRatTraining6moQCRep

### About this repository 
This repository provides R notebooks with the post-quantification analyses of the data presented in the main paper for the first 
large-scale multi-omic multi-tissue endurance exercise training study conducted 
in young adult rats by the Molecular Transducers of Physical Activity Consortium 
(MoTrPAC). Find the [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2022.09.21.508770v2).

Each omic data generated for the paper above had a pipeline to generate quantified features from the raw data (e.g., gene level counts for transcriptomics data quantified from raw RNA-seq reads). In this repository we provide notebooks with the initial analyses of these post-quantification data, which are exploratory in nature. These analyses were used for sanity checks, QA procedures, and identifying flagged samples such as extreme PCA outliers. 

**The reports provided in this repository were originally intended for use within the MoTrPAC consortium, and thus some functionalities that fetch data directly from the internal MoTrPAC cloud systems are not expected to work as is. Nevertheless, we decided to provide these reports because: (1) they contain useful code/analyses and images (e.g., the PCA of each tissue in the transcriptomics report), and (2) most analyses can be reporduced by downloading the data from the [MoTrPAC Data Hub](https://motrpac-data.org/), and modifying the code to take the data from the local directories to which the data were downloaded. The relevant code sections and paths in each report (and its Rmd file) typically use ADDPATH as a token that replaced internal MoTrPAC cloud paths.**

### About MoTrPAC
MoTrPAC is a national research consortium designed to discover and perform 
preliminary characterization of the range of molecular transducers (the 
"molecular map") that underlie the effects of physical activity in humans. 
The program's goal is to study the molecular changes that occur during and after 
exercise and ultimately to advance the understanding of how physical activity 
improves and preserves health. The six-year program is the largest targeted NIH 
investment of funds into the mechanisms of how physical activity improves health 
and prevents disease. See [motrpac.org](https://www.motrpac.org/) and 
[motrpac-data.org](https://motrpac-data.org/) for more details. 

### Repository structure

This repo is organized as follows:

1. `tools` - functions, can source each `.R` file into a session
1. `reports` - the main folder with the reports (`.Rmd` notebooks and their html files)
1. `ext_data` - additional data required for reproducing some analyses (e.g. genome.gtf for the ATAC-seq analysis)

