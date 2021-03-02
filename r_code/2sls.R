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
hes_diag_fname = config$hes_diag_fname
withdrawn_fname <- config$withdrawn_fname
exposure_codes <- config$exposure_codes
outcome_codes <- config$outcome_codes
genotype_fname <- config$genotype_fname
case_covars_output <- config$case_covars_output
gwas_pcs_fname <- config$gwas_pcs_fname
gwas_fname <- config$gwas_fname
p_threshold <- config$p_threshold
clump_threshold <- config$clump_threshold

## Stage 0: Get Data and prepocess

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

## Removing Cancer Cases from the controls
icd10_cancer=c(c(paste0("C0", 0:9), paste0("C", 10:14)), paste0("C", 15:26),
               paste0("C", 30:39), paste0("C", 40:41), paste0("C", 43:44),
               paste0("C", 45:49), paste0("C", 50), paste0("C", 51:58),
               paste0("C", 60:63), paste0("C", 64:68), paste0("C", 69:72),
               paste0("C", 73:75), paste0("C", 76:80), paste0("C", 81:96),
               paste0("C", 97), paste0("D0", 0:9), paste0("D", 37:48))
cancer_cases <- get_hes_cases(hes_diag_fname, withdrawn_fname, icd10_cancer)
case_covars <- case_covars[!(case_covars$eid %in% cancer_cases$eid &
                               case_covars$outcome == 0),]

# drop any remaining nas
case_covars <- na.omit(case_covars)

dim(case_covars)

# convert to factors
case_covars$exposure <- as.factor(case_covars$exposure)
case_covars$outcome <- as.factor(case_covars$outcome)
case_covars$sex <- as.factor(case_covars$sex)
case_covars$immuno <- as.factor(case_covars$immuno)

# save 
saveRDS(case_covars, case_covars_output)

## Stage 1: Regress prs against exposure

case_covars <- readRDS(case_covars_output)

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

# nagelkerke somehow
nagelkerke(exposure_reg)


## Stage 2 regress exposure predictions against the outcome

# get the different function strings
pc_cols <- names(readRDS(gwas_pcs_fname))[-1]
outcome_str <- "outcome ~ "
covars_str <-  paste(c("birth_year", "sex", "immuno", pc_cols), collapse='+')
covars_exp_str <- paste(c(covars_str, "exposure_pred"), collapse='+')

# add the predictions to the dataframe
case_covars$exposure_pred <- exposure_pred

# for only the exposure predictions
outcome_reg <- glm(outcome~exposure_pred, data = case_covars, family = binomial(link = "logit"))
summary(outcome_reg)

# for only the covariates
outcome_covars <- glm(as.formula(c(outcome_str, covars_str)), data = case_covars, family = binomial(link = "logit"))
summary(outcome_covars)

# for the covariates + the exposure predictions
outcome_exp_covars <- glm(as.formula(c(outcome_str, covars_exp_str)), data = case_covars, family = binomial(link = "logit"))
summary(outcome_exp_covars)

# anova after adding the exposure predictions
anova(outcome_covars, outcome_exp_covars)
