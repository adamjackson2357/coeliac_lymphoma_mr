# Run the two-stage least squares

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
library(yaml)
library(ROCR)
library(AER)

# read in the config
config <- read_yaml('../configs/main.yml')
case_covars_output <- config$case_covars_output
prs_output <- config$prs_output

## Read in the risk scores and create prs

# join the data
case_covars <- readRDS(case_covars_output)
prs <- readRDS(prs_output)

# join the prs to the covariates dataframe
case_covars$prs <- prs[match(case_covars$eid, names(prs))]
head(case_covars)

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
plot(exposure_roc, main = "coeliac ~ prs ROC")
abline(a=0, b=1)

# get the auc
exposure_auc= performance(exposure_rocr, measure = "auc")
exposure_auc@y.values

# McFadden's R-squared

# Deviance
null_reg <- glm(exposure~1, data = case_covars, family = binomial(link = "logit"))
anova(null_reg, exposure_reg)

## 2SLS using exposure predictions

# convert to factor
case_covars$outcome <- as.factor(case_covars$outcome)

# logistic regression
outcome_reg <- glm(outcome~exposure_pred, data = case_covars, family = binomial(link = "logit"))
summary(outcome_reg)

## 2SLS using ivreg

# for only exposure
head(case_covars)
tsls <- ivreg(outcome~exposure | prs, data = case_covars)
summary(tsls)

# including covariates
sls_covars <- ivreg(outcome~exposure + sex + age + immuno | prs, data = case_covars)
summary(sls_covars)
