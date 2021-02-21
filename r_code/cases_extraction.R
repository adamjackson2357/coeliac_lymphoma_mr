# Code to extract case eids for a given ICD10 code

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# import libraries
library(yaml)
library(tidyverse)
source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

# Define disease codes
exposure_codes <- config$exposure_codes
outcome_codes <- config$outcome_codes

# Get the exposure cases
exposure <- get_cases(config$hes_diag_fname, config$hes_fname, config$withdrawn_fname, exposure_codes)
saveRDS(exposure, config$exposure_output)

# Get the outcome cases
outcome <- get_cases(config$hes_diag_fname, config$hes_fname, config$withdrawn_fname, outcome_codes)
saveRDS(outcome, config$outcome_output)
