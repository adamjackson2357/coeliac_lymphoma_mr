# Data exploration code

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# read in libraries
library(data.table)
library(openxlsx)
library(tidyverse)

# read in the config
config <- read_yaml('../configs/main.yml')

# read in the case eids
cases=read_rds(config$cases_output)
dim(cases)

# read in the dataset with participant information
covariates=read_rds(config$covariates_output)
dim(covariates)

# Biomarker measurements are stored in a different dataset:
biomarkers=data.frame(fread(config$biomarkers_fname))

# Join cases with covariates and biomarkers
case_covariates <- cases %>%
  left_join(covariates, by="eid") %>%
  left_join(biomarkers, by="eid")
dim(case_covariates)

# get gender counts and proportions
case_covariates %>%
  group_by(sex=X31.0.0) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

