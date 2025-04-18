This repository contains scripts for performing Mendelian Randomization (MR) analysis using the TwoSampleMR R package, followed by post-analysis with MR-PRESSO to detect and correct for pleiotropy. The workflow analyzes two GWAS datasets: one for the exposure and one for the outcome.
Overview
The analysis consists of two main steps:

TwoSampleMR Analysis: We use the TwoSampleMR package to perform two-sample Mendelian Randomization, estimating the causal effect of the exposure on the outcome using genetic variants as instrumental variables.
MR-PRESSO Post-Analysis: After the initial MR analysis, we apply MR-PRESSO (Mendelian Randomization Pleiotropy RESidual Sum and Outlier) to identify horizontal pleiotropy, detect outliers, and obtain pleiotropy-corrected estimates.
Dependencies
To run the scripts, ensure you have the following R packages installed:

TwoSampleMR: For performing two-sample MR analysis.
MRPRESSO: For pleiotropy detection and correction.
Additional dependencies (automatically installed with the above packages): dplyr, ggplot2, etc.

install.packages("TwoSampleMR")
install.packages("MRPRESSO")

Data Requirements
The analysis requires two GWAS summary statistics datasets:

Exposure GWAS: A dataset containing summary statistics for the exposure trait (e.g., SNP IDs, effect sizes, standard errors, p-values).
Outcome GWAS: A dataset containing summary statistics for the outcome trait, formatted similarly.
Ensure the datasets are in a compatible format for TwoSampleMR (e.g., tab-delimited text files with columns for SNP, beta, se, pval, etc.). Example file paths in the scripts are placeholders and should be updated to point to your data.

Workflow
Data Preparation:
Load the exposure and outcome GWAS datasets.
Harmonize the data to ensure SNPs are aligned between exposure and outcome (e.g., matching effect alleles).
TwoSampleMR Analysis:
Extract instrumental variables (IVs) from the exposure GWAS.
Perform MR analysis using methods like Inverse Variance Weighted (IVW), MR-Egger, Weighted Median, and Weighted Mode.
Generate scatter plots, forest plots, and funnel plots to visualize results.
MR-PRESSO Post-Analysis:
Run MR-PRESSO to test for horizontal pleiotropy.
Identify and remove outlier SNPs if significant pleiotropy is detected.
Re-run MR analysis with the corrected dataset to obtain pleiotropy-adjusted estimates.
Output results, including pleiotropy p-values and corrected causal effect estimates.
