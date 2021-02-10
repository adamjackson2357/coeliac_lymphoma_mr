# Code to extract case eids for a given ICD10 code

# Clear variables and set the path
# dev.off()
# rm(list=ls())
# path=dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(path)

# import libraries
library(data.table)
library(tidyverse)
library(yaml)
source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

# Define disease codes
codes <- config$icd10

# pull in the covariates data for participants with cancer
covariates <- get_covariates(config$covariates_fname, config$withdrawn_fname, fields = 40006)

# get the cases
cases <- get_cases(covariates, codes)

# save the cases
saveRDS(cases, config$cases_output)
