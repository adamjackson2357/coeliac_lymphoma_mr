# Count excluded participants

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
library(yaml)
library(data.table)
library(dplyr)

source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

covars_output <- config$covars_output
exposure_output <- config$exposure_output
outcome_output <- config$outcome_output
case_covars_output <- config$case_covars_output
hes_diag_fname <- config$hes_diag_fname
hes_fname <- config$hes_fname
withdrawn_fname <- config$withdrawn_fname
exposure_name <- config$exposure_name
outcome_name <- config$outcome_name
gwas_pcs_fname <- config$gwas_pcs_fname
outcome_codes <- config$outcome_codes
exposure_codes <- config$exposure_codes
fields <- config$fields_fname
covars_fname <- config$covars_fname
outcome_codes


hes_outcome<- get_hes_cases(hes_diag_fname, withdrawn_fname, outcome_codes)
head(hes_outcome)

# get the covariate column names
column_ids <- data.frame(fread(covars_fname, nrows=1))

# get the columns ids for the given fields
column_ids <- extract_column_ids(column_ids, fields)

# Extracting required columns from dataset
covariates <- data.frame(fread(covars_fname, select=column_ids))

# lots of "" instead of Nans, replace these
covariates[covariates == ""] <- NA

# first count
dim(covariates)

# remove participants which have withdrawn
withdrawn=read.csv(withdrawn_fname)[,1]
length(withdrawn)
all_participants <- subset(covariates, !(eid %in% withdrawn))
dim(all_participants)

# get the case_covariates
case_covars <- get_case_covars(covars_fname, hes_diag_fname, withdrawn_fname,
                            fields_fname, exposure_codes, outcome_codes)
dim(case_covars)
head(case_covars)

table(case_covars$outcome)

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

table(case_covars$outcome)

head(case_covars)
sapply(case_covars, function(y) sum(length(which(is.na(y)))))
case_covars <- na.omit(case_covars)
table(case_covars$outcome)


case_covars <- readRDS(case_covars_output)
table(case_covars$outcome)
