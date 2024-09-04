# Meta-Analysis of Epigenome-Wide Association Study (EWAS) on Ultra-Processed Food Consumption and DNA Methylation in European Children
## Description
This project includes scripts used to perform various steps in an EWAS, focusing on the association between ultra-processed food (UPF) consumption and DNA methylation (DNAm) in European children. The following is a summary of each step:

### 1 - Descriptive and Association Analysis:
This script generates descriptive statistics and association analyses with cell types across all cohorts, providing tables and plots. It is run prior to the EWAS to examine potential associations with cell type composition.

### 2 - Run EWAS:
Implements a robust linear model (RLM) to assess the association between UPF and DNAm. The output consists of winsorized beta values, with UPF as the primary exposure variable.

### 3 - Quality Control:
This step reads the EWAS results and performs the following:

### 4 - Computes adjusted p-values (if not already done).
Prepares the data for quality control using the EASIER package.
Runs quality control using the EASIER package to remove unreliable probes and ensure data quality.
Meta-Analysis of DMPs:
This script runs a meta-analysis of differentially methylated positions (DMPs) across cohorts using the GWAMA software.

### 5 - Analyze Meta DMPs:
Checks the quality of the meta-analysis by visualizing QQ plots, computing lambda values, and selecting the top CpGs. It also generates Manhattan and forest plots for selected CpGs.

### 6 - Correlation Matrix for DMRs:
Calculates correlations between each CpG site and the subsequent 20 CpGs in the genome. This information is anonymized and safe to share with external users as it does not reveal individual-level data.

### 7 - Meta-Analysis of DMRs:
Performs a meta-analysis by combining the correlation matrices of differentially methylated regions (DMRs) across cohorts after running the EWAS for DMPs.

### 8 - Enrichment of DMPs:
Conducts functional enrichment analysis of selected CpGs based on a specified threshold to gain biological insight.


## R Session Information

R version: 4.3.0 (2023-04-21)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Rocky Linux 8.7 (Green Obsidian)

### Key Packages Used:
msigdbr: 7.5.1
IlluminaHumanMethylationEPICanno: 0.6.0
IlluminaHumanMethylation450kanno.ilmn12.hg19: 0.6.1
minfi: 1.44.0
EASIER: 0.1.2.28
missMethyl: 1.32.1
tidyverse: 2.0.0
brgeEnrich: 0.1.0
GWAMA: Required for meta-analysis. It has to be installed in the cluser where the meta-analysis run. 
ggplot2: 3.4.2
GenomicRanges: 1.50.2
dmrff: 1.1.1 