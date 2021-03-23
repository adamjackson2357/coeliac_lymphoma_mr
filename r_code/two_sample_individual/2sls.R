# Run the two-stage least squares

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
  library(ggpubr)
  library(forestmodel)
  library(xtable)
  library(patchwork)
  library(flextable)
  library(caret)
})
source("../utils/extraction.R")
source("../2sample_individual/prs.R")
source("../2sample_individual/model_processing.R")

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

case_covars <- get_features(covars_fname, hes_diag_fname, withdrawn_fname,
                         fields_fname, exposure_codes, outcome_codes,
                         genotype_fname, gwas_fname, clump_threshold,
                         gwas_pcs_fname)
saveRDS(case_covars, case_covars_output)
case_covars <- readRDS(case_covars_output)


## Stage 1: Regress prs against exposure
case_covars <- readRDS(case_covars_output)
case_covars <- gluten_free(case_covars, covars_fname, withdrawn_fname)
head(case_covars)
str(case_covars)

# Coeliac logistic regression
exposure_reg <- glm(exposure~prs, data = case_covars, family = binomial(link = "logit"))
exposure_pred <- predict(exposure_reg, type="response", newdata=case_covars)
exposure_rocr <- prediction(exposure_pred, case_covars$exposure, label.ordering = NULL)
performance(exposure_rocr, measure = "auc")@y.values
exposure_roc <- performance(exposure_rocr, measure="tpr", x.measure="fpr")
exposure_values = data.frame(fpr=exposure_roc@x.values[[1]],
                             fnr=exposure_roc@y.values[[1]])
exposure_values$Exposure <- "Coeliac"


# Gluten free logistic regression
gluten_reg <- glm(Gluten_free~prs, data = case_covars, family = binomial(link = "logit"))
gluten_pred <- predict(gluten_reg, case_covars, type = "response")
gluten_rocr <- prediction(gluten_pred, case_covars$Gluten_free, label.ordering = NULL)
performance(gluten_rocr, measure = "auc")@y.values
gluten_roc <- performance(gluten_rocr, measure="tpr", x.measure="fpr")
gluten_values = data.frame(fpr=gluten_roc@x.values[[1]],
                           fnr=gluten_roc@y.values[[1]])
gluten_values$Exposure <- "Gluten Free"

# get the prs roc
fpr_fnr <- rbind(exposure_values, gluten_values)
prs_roc <-  fpr_fnr %>%
  ggplot(aes(x = fpr, y = fnr, colour=Exposure)) +
  geom_line() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "black") +
  theme_minimal() +
  xlim(0, 1) + ylim(0, 1) +
  xlab("False Positive Rate") + ylab("False Negative Rate")
prs_roc
ggsave("../../figures/prs_roc.png", prs_roc, width=10, height=5)

## Stage 2 regress exposure predictions against the outcome

case_covars <- readRDS(case_covars_output)
case_covars <- feature_recoding(case_covars)
case_covars$Coeliac_Prediction <- exposure_pred
head(case_covars)

# get the different function strings
pc_cols <- names(readRDS(gwas_pcs_fname))[2:11]
outcome_str <- "NHL ~ "
covars_str <-  paste(c("Sex", "Age", pc_cols), collapse='+')
covars_exp_str <- paste(c(covars_str, "Coeliac_Prediction"), collapse='+')

exp_model <- glm(as.formula(c(outcome_str, covars_exp_str)), data=case_covars, family=binomial(link="logit"))
summary(exp_model)
exp(exp_model$coefficients)
exp(confint(exp_model))
covars_list = c("Sex", "Age", "Coeliac_Prediction")
format_list <- forest_model_format_options(banded = FALSE)
sls_forest <- forest_model(exp_model, factor_separate_line=TRUE,
                            covariates=covars_list, format_options = format_list)
ggsave(plot = sls_forest, filename="../../figures/sls_forest.png", width=10, height=5)

# for only the covars predictions
covars_model <- glm(as.formula(c(outcome_str, covars_str)), data=case_covars, family=binomial(link="logit"))
summary(covars_model)
anova(covars_model, exp_model)

