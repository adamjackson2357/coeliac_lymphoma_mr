# Code to extract case eids for a given ICD10 code

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# import libraries
library(yaml)
library(tidyverse)
source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

# Define disease codes
codes <- config$risk_factors

# get the cases
cases <- get_cases(config$hes_diag_fname, config$hes_fname, config$withdrawn_fname, codes)

# save the cases
saveRDS(cases, config$risk_factor_output)
