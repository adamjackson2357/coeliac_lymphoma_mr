# Extract the covariates data

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(data.table)
source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

# get the fields
fields <- unname(unlist(read.table(config$fields_fname, header=FALSE)))

# get the covariates
covariates <- get_covariates(config$covariates_fname,
                             config$withdrawn_fname,
                             fields)

# save
saveRDS(covariates, config$covariates_output)
