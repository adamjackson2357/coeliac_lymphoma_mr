# TDS Project

## Introduction

Group 8 analysis for TDS Module.

Using Mendelian Randomisation to assess the causative effects of coeliac disease on
non-hodgkins lymphoma.

## Quick start

1. Clone the project from [https://github.com/adamjackson2357/coeliac_lymphoma_mr]

## Prerequisites

r_packages:
- devtools
- yaml
- data.table
- tidyverse
- tableone
- forestmodel
- ieugwasr
- EBPRS
- ROCR
- ggplot2
- ggpubr
- xtable
- patchwork
- flextable
- DiagrammeR

The following UK biobank files placed in the data folder
- "./data/ukb26390.csv"
- "./data/w19266_20200204.csv"
- "./data/hesin_diag.txt"
- "./data/hesin.txt"
- "./data/genetic_data_extracted.rds"
- "./data/GWAS_PCs.rds"

## Documentation

The documentation for this project is this README.

## Analysis

Observational Analysis
- Codes to create the descriptive statistics and measure association between coeliac disease and NHL

Two Sample Individual MR
- Codes to create the PRS and then run the two-stage least squares analysis

Two Sample Summary MR
- Codes to run the two sample summary MR analysis

GWAS
- Codes to run the GWAS

## Support team

- Ahmed Abdulaal
- Aimee Tham
- Irene Sugurova
- Adam Jackson
