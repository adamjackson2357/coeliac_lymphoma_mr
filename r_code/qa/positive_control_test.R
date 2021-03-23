# PRS Positive Control Test

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
})
source("extraction.R")
source("prs.R")

# read in the config
config <- read_yaml('../configs/main.yml')
covars_fname <- config$covars_fname
withdrawn_fname <- config$withdrawn_fname
case_covars_output <- config$case_covars_output
genotype_fname <- config$genotype_fname
gwas_fname <- config$gwas_fname
p_threshold <- config$p_threshold
clump_threshold <- config$clump_threshold

diet <- get_covariates(covars_fname, withdrawn_fname, "20086")
gluten_free <- diet[,c(1, grep(pattern = ".0$", colnames(diet), ))]
gluten_free <- subset(gluten_free, X20086.0.0 == 8)
case_covars$exposure <- ifelse(case_covars$eid %in% gluten_free$eid, 1, 0)

## Get the Data

# read in the case covariate data
case_covars <- readRDS(case_covars_output)

# get gluten free
diet <- get_covariates(covars_fname, withdrawn_fname, "20086")

# filter for people that have been gluten free at some point
gluten_free <- diet[,c(1, grep(pattern = ".0$", colnames(diet), ))]
gluten_free <- subset(gluten_free, X20086.0.0 == 8)
head(gluten_free)
# gluten_free <- gluten_free %>%
#   filter_all(any_vars(. %in% 8))

case_covars$exposure <- ifelse(case_covars$eid %in% gluten_free$eid, 1, 0)
head(case_covars)

prs <- get_prs(genotype_fname, gwas_fname, 5*10^-8, clump_threshold)
case_covars$prs <- prs[match(case_covars$eid, names(prs))]

## Predict the PRS on being gluten free

table(case_covars$exposure)

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
