# SNP Data Processing and PCA Analysis

This repository provides a reproducible workflow for processing SNP genotype data using PLINK and visualizing population structure through Principal Component Analysis (PCA) in R.

## Overview

The script automates the following steps:
- Converts raw VCF files to PLINK formats (PED/MAP and binary).
- Performs linkage disequilibrium (LD)-based SNP pruning to reduce marker redundancy.
- Computes a genetic distance matrix from the pruned dataset.
- Performs PCA using classical multidimensional scaling.
- Generates PCA plots with and without k-means clustering to explore genetic structure and subpopulation patterns.

This workflow is useful for population structure analysis, diversity assessment, and as a preprocessing step for GWAS or genomic prediction studies.
