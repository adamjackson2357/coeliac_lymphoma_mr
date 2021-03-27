# Run the two-stage least squares

# Clear variables and set the path
# dev.off()
# rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(ROCR)
  library(ggplot2)
  library(data.table)
  library(ggpubr)
})
source("../utils/extraction.R")
source("../two_sample_individual/prs.R")
source("../two_sample_individual/model_processing.R")
source("../two_sample_individual/create_roc.R")

# read in the config
config <- read_yaml('../../configs/main.yml')
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
clump_threshold <- config$clump_threshold

## Stage 0: Get Data and preprocess

# case_covars <- get_features(covars_fname, hes_diag_fname, withdrawn_fname,
#                             fields_fname, exposure_codes, outcome_codes,
#                             genotype_fname, gwas_fname, clump_threshold,
#                             gwas_pcs_fname)
# saveRDS(case_covars, case_covars_output)
# case_covars <- readRDS(case_covars_output)

## Stage 1: Regress prs against exposure
case_covars <- readRDS(case_covars_output)
case_covars <- gluten_free(case_covars, covars_fname, withdrawn_fname)

# Coeliac logistic regression
exposure_reg <- glm(exposure~prs, data = case_covars, family = binomial(link = "logit"))
exposure_values <- get_fpr_fnr(case_covars, exposure_reg, 'exposure', 'Coeliac')

# Gluten free logistic regression
gluten_reg <- glm(Gluten_free~prs, data = case_covars, family = binomial(link = "logit"))
gluten_values <- get_fpr_fnr(case_covars, gluten_reg, 'Gluten_free', "Gluten Free")

# get the exposure prs roc
exposure_fpr_fnr <- rbind(exposure_values, gluten_values)
rm(exposure_values)
rm(gluten_values)
exposure_roc <- create_roc_curve(exposure_fpr_fnr)
rm(exposure_fpr_fnr)
ggsave("../../figures/exposure_roc.png", exposure_roc, width=10, height=5)

## Stage 2 regress exposure predictions against the outcome

# get the different function strings
pc_cols <- names(readRDS(gwas_pcs_fname))[2:11]
outcome_str <- "outcome ~ "
covars_str <-  paste(c("sex", "age", pc_cols), collapse='+')
covars_prs_str <- paste(c(covars_str, "prs"), collapse='+')

# NHL logistic regression
prs_reg <- glm(outcome~prs, data = case_covars, family = binomial(link = "logit"))
prs_values <- get_fpr_fnr(case_covars, prs_reg, 'outcome', 'PRS')

# Covars logistic regression
covars_reg <- glm(as.formula(c(outcome_str, covars_str)), data = case_covars, family = binomial(link = "logit"))
covars_values <- get_fpr_fnr(case_covars, covars_reg, 'outcome', 'Covariates')

# Covars + NHL logistic regression
covars_prs_reg <- glm(as.formula(c(outcome_str, covars_prs_str)), data = case_covars, family = binomial(link = "logit"))
covars_prs_values <- get_fpr_fnr(case_covars, covars_prs_reg, 'outcome', 'Covarariates + PRS')

# get the outcome prs roc
outcome_fpr_fnr <- rbind(covars_prs_values, covars_values, prs_values)
rm(prs_values)
rm(covars_values)
rm(covars_prs_values)
outcome_roc <- create_roc_curve(outcome_fpr_fnr)
rm(outcome_fpr_fnr)
ggsave("../../figures/outcome_roc.png", outcome_roc, width=10, height=5)

# Print both side by side
rocs <- ggarrange(
  exposure_roc, outcome_roc, ncol=2, nrow=1,
  common.legend = FALSE, legend = "bottom"
)
ggsave(plot = rocs, filename="../../figures/rocs.png")
