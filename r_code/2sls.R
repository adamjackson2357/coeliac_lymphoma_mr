# Run the two-stage least squares

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
library(yaml)
library(dplyr)
library(ieugwasr)
library(EBPRS)
library(ROCR)
source("extraction.R")
source("prs.R")

# read in the config
config <- read_yaml('../configs/main.yml')
fields_fname <- config$fields_fname
covars_fname <- config$covars_fname
withdrawn_fname <- config$withdrawn_fname
covars_output <- config$covars_output
hes_diag_fname = config$hes_diag_fname
withdrawn_fname <- config$withdrawn_fname
exposure_codes <- config$exposure_codes
outcome_codes <- config$outcome_codes
exposure_output <- config$exposure_output
outcome_output <- config$outcome_output
genotype_fname <- config$genotype_fname
ivs_output <- config$ivs_output
case_covars_output <- config$case_covars_output
gwas_pcs_fname <- config$gwas_pcs_fname
gwas_fname <- config$gwas_fname
p_threshold <- config$p_threshold
clump_threshold <- config$clump_threshold

## Case-covars

case_covars <- get_case_covars(covars_fname, hes_diag_fname, withdrawn_fname,
                            fields_fname, exposure_codes, outcome_codes)

# recoding
old_names <- c("X34.0.0", "X31.0.0")
new_names <- c("birth_year", "sex")
case_covars <- case_covars %>% rename_at(vars(old_names), ~ new_names)

# get participants that are immuno-compromised
immuno <- data.frame(fread(hes_diag_fname, select = c("eid", "diag_icd10")))
immuno <- subset(immuno, startsWith(immuno$diag_icd10, "D8"))
case_covars$immuno <- ifelse(case_covars$eid %in% immuno$eid, 1, 0)

# get the prs and join to the case covariates dataframe
prs <- get_prs(genotype_fname, gwas_fname, 5*10^-8, clump_threshold)
case_covars$prs <- prs[match(case_covars$eid, names(prs))]

# get the gwas pcs and join
gwas_pcs <- readRDS(gwas_pcs_fname)
case_covars <- inner_join(case_covars, gwas_pcs, by=c("eid"))

# save 
saveRDS(case_covars, case_covars_output)

## Stage 0: Get data and preprocess

case_covars <- readRDS(case_covars_output)
pc_cols <- names(readRDS(gwas_pcs_fname))[-1]

case_covars$exposure <- as.factor(case_covars$exposure)
case_covars$outcome <- as.factor(case_covars$outcome)
case_covars$sex <- as.factor(case_covars$sex)
case_covars$immuno <- as.factor(case_covars$immuno)

## Stage 1: Regress prs against exposure

# logistic regression
exposure_reg <- glm(exposure~prs, data = case_covars, family = binomial(link = "logit"))
summary(exposure_reg)

# predict on exposure
exposure_pred <- predict(exposure_reg, case_covars, type = "response")

# create the ROCR object
exposure_rocr <- prediction(exposure_pred, case_covars$exposure, label.ordering = NULL)

# ROC Curve
exposure_roc <- performance(exposure_rocr, measure="tpr", x.measure="fpr")
plot(exposure_roc, main = "coeliac ~ prs ROC")
abline(a=0, b=1)

# get the auc
exposure_auc= performance(exposure_rocr, measure = "auc")
exposure_auc@y.values

# Pseudo R-squared
1 - exposure_reg$deviance / exposure_reg$null.deviance

# Deviance
null_reg <- glm(exposure~1, data = case_covars, family = binomial(link = "logit"))
anova(null_reg, exposure_reg)

## Stage 2 regress exposure predictions against the outcome

# for only the exposure predictions
outcome_reg <- glm(outcome~exposure_pred, data = case_covars, family = binomial(link = "logit"))
summary(outcome_reg)

# for only the covariates
outcome_covars <- glm(outcome~sex+age+immuno, data = case_covars, family = binomial(link = "logit"))
summary(outcome_covars)

# for the covariates + the exposure predictions
outcome_exp_covars <- glm(outcome~sex+age+immuno+exposure_pred, data = case_covars, family = binomial(link = "logit"))
summary(outcome_exp_covars)

# anova after adding the exposure predictions
anova(outcome_covars, outcome_exp_covars)

