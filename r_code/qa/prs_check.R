# Harmonisation code creation

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(ieugwasr)
  library(EBPRS)
  library(ROCR)
  library(ggplot2)
  library(forestmodel)
})
source("extraction.R")
source("prs.R")

# read in the config
config <- read_yaml('../configs/main.yml')
fields_fname <- config$fields_fname
covars_fname <- config$covars_fname
withdrawn_fname <- config$withdrawn_fname
hes_diag_fname = config$hes_diag_fname
withdrawn_fname <- config$withdrawn_fname
exposure_codes <- config$exposure_codes
outcome_codes <- config$outcome_codes
genotype_fname <- config$genotype_fname
case_covars_output <- config$case_covars_output
gwas_pcs_fname <- config$gwas_pcs_fname
gwas_fname <- config$gwas_fname
p_threshold <- 5*10^-8
clump_threshold <- config$clump_threshold
allele_harmonisation_fname <- config$allele_harmonisation_fname

allele_harmonisation_fname


# read in the data
geno <- readRDS(genotype_fname)
ivs <- get_ivs(gwas_fname, p_threshold, clump_threshold)
harmon <- read.csv(allele_harmonisation_fname)[, c("SNP", "V2")]
harmon
ivs <- ivs %>%
  inner_join(harmon, by=c("rsid"="SNP"))
snps <- ivs$rsid
log_or <- log(ivs$OR)
harmon_alleles <- ivs$V2

# preprocess
geno <- preprocessing(geno, snps)

# harmonise the log odds ratios
log_or <- harmonise(log_or, minor_alleles)

# get the prs score
prs <- mean_prs(geno, log_or)

