# Run the two-stage least squares

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
library(yaml)
library(ROCR)

# read in the config
config <- read_yaml('../configs/main.yml')
prs_output <- config$prs_output
case_covars_output <- config$case_covars_output

## Read in the PRS and join with covariates

# read in the PRS Data
prs <- readRDS(prs_output)

# join the case-covariates dataframe
case_covars <- readRDS(case_covars_output)

# join the prs to the covariates dataframe
case_covars$prs <- prs[match(case_covars$eid, names(prs))]

# check number of missing values
length(case_covars$prs)
length(na.omit(case_covars$prs))

# reduce down the dataframe
case_covars <- case_covars[, c("eid", "sex", "exposure", "outcome", "age", "immuno", "prs")]
case_covars <- na.omit(case_covars)
length(na.omit(case_covars$prs))

## Regress against exposures disease

# convert to factor
case_covars$exposure <- as.factor(case_covars$exposure)

# logistic regression
exposure_reg <- glm(exposure~prs, data = case_covars, family = binomial(link = "logit"))
summary(exposure_reg)

# predict on exposure
exposure_pred <- predict(exposure_reg, case_covars, type = "response")

# create the ROCR object
exposure_rocr <- prediction(exposure_pred, case_covars$exposure, label.ordering = NULL)

# ROC Curve
exposure_roc <- performance(exposure_rocr, measure="tpr", x.measure="fpr")
plot(exposure_roc)
abline(a=0, b=1)

# get the auc
exposure_auc= performance(exposure_rocr, measure = "auc")
exposure_auc@y.values

## Regress the exposure predictions against lymphoma

# convert to factor
case_covars$outcome <- as.factor(case_covars$outcome)

# logistic regression
outcome_reg <- glm(outcome~exposure_pred, data = case_covars, family = binomial(link = "logit"))
summary(outcome_reg)

# predict on exposure
outcome_pred <- predict(outcome_reg, case_covars, type = "response")

# create the ROCR object
outcome_rocr <- prediction(outcome_pred, case_covars$outcome, label.ordering = NULL)

# ROC Curve
outcome_roc <- performance(outcome_rocr, measure="tpr", x.measure="fpr")
plot(outcome_roc)
abline(a=0, b=1)

# get the auc
outcome_auc = performance(outcome_rocr, measure = "auc")
outcome_auc@y.values
