# RNA-seq Analysis for Mtmr6 Knockout Mouse Model

This repository contains the custom R scripts and processed data used for the transcriptomic analysis in our manuscript investigating the role of *Mtmr6* in male fertility maintenance.

## System Requirements
* **Software**: R version 4.4.3 or later.
* **Operating System**: Windows, macOS, or Linux.

## Required R Packages
Ensure the following standard open-source packages are installed before running the scripts:
* `DESeq2`
* `clusterProfiler`
* `org.Mm.eg.db`
* `ggplot2`
* `pheatmap`
* `enrichplot`
* `eulerr`

## Repository Structure
* `data/`: Contains processed count matrices and differential expression results (.csv files) for heterozygous (Het) and knockout (KO) mouse models.
* `scripts/`: Contains the R scripts for data analysis and visualization.
    * `main_analysis.R`: Primary script for differential gene expression analysis and Gene Ontology (GO) enrichment.
    * `volcano_plots_analysis.R`: Script for generating volcano plots to visualize specific up- and down-regulated genes (e.g., *Espn, Pea15a, Bik, Defb36*).

## Instructions for Use
1. Download or clone this repository to your local machine.
2. Set the working directory in RStudio to the root folder of this repository.
3. Execute the scripts in the `scripts/` directory to reproduce the transcriptomic analysis and figures presented in the manuscript.
