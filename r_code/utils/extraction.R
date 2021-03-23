## Set of functions for extracting data

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
#' @param covars_fname filename for all the covariates
#' @param withdrawn_fname filename for withdrawn participants
#' @param fields list of the column ids
#' @return Covariates dataframe
get_covariates <- function(covars_fname, withdrawn_fname, fields) {
  
  # get the covariate column names
  column_ids <- data.frame(fread(covars_fname, nrows=1))
  
  # get the columns ids for the given fields
  column_ids <- extract_column_ids(column_ids, fields)
  
  # Extracting required columns from dataset
  covariates <- data.frame(fread(covars_fname, select=column_ids))
  
  # lots of "" instead of Nans, replace these
  covariates[covariates == ""] <- NA
  
  # remove participants which have withdrawn
  withdrawn=read.csv(withdrawn_fname)[,1]
  covariates <- subset(covariates, !(eid %in% withdrawn))
  
  return(covariates)
}


#' Extract the cases from the hes dataset
#' 
#' Filters the hes dataframe for a set of icd10 codes
#' 
#' @param hes_diag_fname hes fname
#' @param withdrawn_fname filename for withdrawn participants
#' @param codes list of icd10 disease codes
#' @return eids for a given set of icd10 codes
get_hes_cases <- function(hes_diag_fname, withdrawn_fname, codes){
  
  # get the hes_diag data
  hes_diag <- data.frame(fread(hes_diag_fname,
                               select = c("eid", "diag_icd10")))
  
  # filter on icd 10 code and return distinct eids
  cases <- hes_diag %>%
    filter(diag_icd10 %in% codes) %>%
    select(eid) %>%
    distinct()
  
  # remove participants which have withdrawn
  withdrawn=read.csv(withdrawn_fname)[,1]
  cases <- subset(cases, !(eid %in% withdrawn))
  
  return(cases)
}


#' Extract the cases from the covariates dataset
#' 
#' Select specific columns from the covariates dataset and filters for a set
#' of icd10 codes
#' 
#' @param covars_fname covariates fname
#' @param withdrawn_fname filename for withdrawn participants
#' @param fields list of the column ids
#' @param codes list of icd10 disease codes
#' @return eids for a given set of icd10 codes
get_covars_cases <- function(covars_fname, withdrawn_fname, fields, codes){
  
  column_ids <- data.frame(fread(covars_fname, nrows=1))
  
  # select only the cancer column
  covars <- get_covariates(covars_fname, withdrawn_fname, fields)
  
  # select anyone that has one of the icd10 codes
  cases <- covars %>%
    filter_all(any_vars(. %in% codes)) %>%
    select(eid) %>%
    distinct()
  
  # remove participants which have withdrawn
  withdrawn=read.csv(withdrawn_fname)[,1]
  cases <- subset(cases, !(eid %in% withdrawn))
  
  return(cases)
}


get_case_covars <- function(covars_fname, hes_diag_fname, withdrawn_fname,
                            fields_fname, exposure_codes, outcome_codes){
  
  # get the fields
  fields <- unname(unlist(read.table(config$fields_fname, header=FALSE)))
  
  # get the covariates
  covars <- get_covariates(covars_fname,
                           withdrawn_fname,
                           fields)
  
  # Get the exposure cases
  exposure <- get_hes_cases(hes_diag_fname, withdrawn_fname, exposure_codes)
  exposure$exposure <- 1
  
  # Get the outcome cases
  hes_outcome<- get_hes_cases(hes_diag_fname, withdrawn_fname, outcome_codes)
  covars_outcome <- get_covars_cases(covars_fname, withdrawn_fname, c("40006"), outcome_codes)
  outcome <- full_join(hes_outcome, covars_outcome, by=c("eid"))
  outcome$outcome <- 1
  head(outcome)
  
  # join the three dataframes to create case-covars
  case_covars <- covars %>%
    left_join(exposure, by="eid") %>%
    left_join(outcome, by="eid", suffix=c("_exposure", "_outcome")) %>%
    mutate(exposure = if_else(is.na(exposure), 0, exposure)) %>%
    mutate(outcome = if_else(is.na(outcome), 0, outcome))
  
  return(case_covars)
}


#' Extract the cases from the hes datasets including the episode durations
#' 
#' Filters the covariates dataframe for a set of icd10 codes and melts it to
#' a long dataframe. Then extracts the time points as a new columns
#' 
#' @param hes_diag_fname hes fname
#' @param hes_fname hes fname with date of diagnosis
#' @param withdrawn_fname filename for withdrawn participants
#' @param codes list of icd10 disease codes
#' @return Covariates dataframe
get_cases_duration <- function(hes_diag_fname, hes_fname, withdrawn_fname, codes){
  
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
