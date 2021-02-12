## Set of functions for extracting data

library(data.table)
library(tidyverse)

#' Extract the column ids from a list of field ids
#' 
#' @param data Dataframe containing all the field_ids, e.g. ukb26390.csv
#' @param fields list of field IDs that we want to extract
#' @return List of the column ids
extract_column_ids <- function(data, fields) {
  
  # Find the position of the column called eid and save it 
  column_ids=grep("eid", colnames(data))
  
  # Create an empty vector 
  found_fieldids=NULL
  
  # Loop through the field IDs and return positions of the columns of interest
  for (k in 1:length(fields)){
    mygrep=grep(paste0("X",fields[k],"."), fixed=TRUE, colnames(data))
    if (length(mygrep)>0){
      found_fieldids=c(found_fieldids, fields[k])
    }
    column_ids=c(column_ids, mygrep)
  }
  
  return(column_ids)
}


#' Extract the covariates dataframe for a set of field ids
#' 
#' Gets the column names, filters by field id, replaces empty strings
#' and removes those that have withdrawn
#' 
#' @param covariates_fname filename for all the covariates
#' @param withdrawn_fname filename for withdrawn participants
#' @param fields list of the column ids
#' @return Covariates dataframe
get_covariates <- function(covariates_fname, withdrawn_fname, fields) {
  
  # get the covariate column names
  column_ids <- data.frame(fread(covariates_fname, nrows=1))
  
  # get the columns ids for cancer only
  column_ids <- extract_column_ids(column_ids, fields)
  
  # Extracting required columns from dataset
  covariates <- data.frame(fread(covariates_fname, select=column_ids))
  
  # lots of "" instead of Nans, replace these
  covariates[covariates == ""] <- NA
  
  # remove participants which have withdrawn
  withdrawn=read.csv(withdrawn_fname)[,1]
  covariates <- subset(covariates, !(eid %in% withdrawn))
  
  return(covariates)
}


#' Extract the cases from the hes datasets
#' 
#' Filters the covariates dataframe for a set of icd10 codes and melts it to
#' a long dataframe. Then extracts the time points as a new columns
#' 
#' @param hes_diag_fname hes fname
#' @param hes_fname hes fname with date of diagnosis
#' @param withdrawn_fname filename for withdrawn participants
#' @param codes list of icd10 disease codes
#' @return Covariates dataframe
get_cases <- function(hes_diag_fname, hes_fname, withdrawn_fname, codes){
  
  # get the hes_diag data
  hes_diag <- data.frame(fread(hes_diag_fname,
                               select = c("eid", "ins_index", "diag_icd10")))
  
  # filter on cases
  hes_diag <- subset(hes_diag, diag_icd10 %in% codes)
  
  # remove participants which have withdrawn
  withdrawn=read.csv(withdrawn_fname)[,1]
  hes_diag <- subset(hes_diag, !(eid %in% withdrawn))
  
  # get the hes data
  hes <- data.frame(fread(hes_fname,
                          select = c("eid", "ins_index", "epistart", "epiend")))
  
  # inner join
  cases <- hes_diag %>%
    inner_join(hes, by = c("eid", "ins_index"))
  
  return(cases)
}
