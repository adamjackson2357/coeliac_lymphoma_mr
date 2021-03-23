# Check the number of lymphoma cases

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(yaml)
library(dplyr)
library(data.table)

source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

hes_diag_fname <- config$hes_diag_fname
outcome_codes <- config$outcome_codes
withdrawn_fname <- config$withdrawn_fname
hes_fname <- config$hes_fname
outcome_output <- config$outcome_output
covars_fname <- config$covars_fname


#### METHOD 1: USING HES

# get the hes_diag data
hes_diag <- data.frame(fread(hes_diag_fname,
                             select = c("eid", "ins_index", "diag_icd10")))

# filter on cases
hes_diag <- subset(hes_diag, diag_icd10 %in% outcome_codes)

# remove participants which have withdrawn
withdrawn=read.csv(withdrawn_fname)[,1]
hes_diag <- subset(hes_diag, !(eid %in% withdrawn))
length(hes_diag$eid)
length(unique(hes_diag$eid))

#### METHOD 2: USING COVARIATES DATAFRAME


# get the biobank data
bb=data.frame(fread("../data/ukb26390.csv"))

# remove participants which have withdrawn
withdrawn=read.csv(config$withdrawn_fname)[,1]
bb <- subset(bb, !(eid %in% withdrawn))


colnames(data_frame)[grep(field_id,colnames(data_frame))]

cancer_cols <- extract_column_ids(bb, "40006")

bb <- data.frame(fread("../data/ukb26390.csv", select=cancer_cols))
head(bb)


for (i in 1:length(colnames(data_frame)[grep(field_id,colnames(data_frame))])){
  print()
}


# get unique eids back for a given disease code
get_biobank_cases <- function(data_frame, field_id, icd_codes) {
  
  #Create empty dataframe with same column names as biobank data in the global environment
  cases <<- data_frame[FALSE,]
  
  #Loop through all time points for which ICD10 codes were collected 
  for (i in 1:length(colnames(data_frame)[grep(field_id,colnames(data_frame))])) {
    #Create variable holding current column name
    column <- paste0("X", field_id,".", i-1, ".0")
    
    #Find everyone with the specified ICD10 codes at the current time point
    disease_at_time_i <- filter(data_frame, !!as.symbol(column) %in% icd_codes)
    
    #Add as rows to the cases dataframe
    cases <<- rbind(cases, disease_at_time_i)
  }
  
  #Only keep unique participants
  cases <<- cases %>%
    distinct(eid)
  #If you want to keep entire datafram = distinct(eid, .keep_all = TRUE)
  
  #Adding column with disease status 
  cases["NHL"] <<- 1
  
  #Create controls dataframe (everyone without the disease(s))
  controls <<- subset(data_frame, !(eid %in% cases[,1]))
  controls <<- controls %>%
    distinct(eid)
  
  #Adding column with disease status 
  controls["NHL"] <<- 0
  
  #Concatenate cases and controls into new dataframe with just eid and disease status
  case_control <<- rbind(cases, controls)
}

#Calling function
get_biobank_cases(bb, "40006", icd10)

## METHOD 3: QUICKER VERSION OF 2
covars_fname

column_ids <- data.frame(fread(covars_fname, nrows=1))

# get the covariates
covars <- get_covariates(covars_fname, withdrawn_fname, c("40006"))

# select anyone that has one of the lymphoma codes
outcome <- covars %>% filter_all(any_vars(. %in% outcome_codes))

length(outcome$eid)

## Compare the overlap between the methods

# eids from the covariates
covars <- readRDS("/rds/general/project/hda_students_data/live/Group8/Ahmed/General/data/cases_controls.rds")
covars <- subset(covars, NHL == 1)
covars_eids <- covars %>%
  select(eid) %>%
  mutate(covars = 1)
head(covars_eids)

# eids from hes data
hes <- readRDS(outcome_output)
hes_eids <- hes %>%
  select(eid) %>%
  distinct() %>%
  mutate(hes = 1)

# join
all_eids <- covars_eids %>%
  full_join(hes_eids, by=c("eid"))

# fill nas
all_eids[is.na(all_eids)] <- 0

apply(all_eids, MARGIN = 2, FUN = function(x) sum(length(which(is.na(x)))))
addmargins(table(covars = all_eids$covars, hes = all_eids$hes))


## Checking the time of first diagnosis code

# get the time of first diagnosis
get_first_diag <- function(df) {
  df <- df %>%
    mutate(epistart = as.Date(epistart, format = "%d/%m/%Y")) %>%
    arrange(eid, epistart) %>%
    group_by(eid) %>%
    mutate(diag = first(epistart)) %>%
    ungroup() %>%
    filter(epistart == diag |is.na(epistart)) %>%
    select(eid, diag_icd10, diag)
  return(df)
}

# read in the outcome data
outcome <- readRDS(outcome_output)

pre_eids <- outcome %>%
  select(eid) %>%
  distinct() %>%
  mutate(pre = 1)

head(pre_eids)
dim(pre_eids)

# run first diagnosis code
outcome <- get_first_diag(outcome)

post_eids <- outcome %>%
  select(eid) %>%
  distinct() %>%
  mutate(post = 1)

head(outcome)
dim(post_eids)

all_eids <- pre_eids %>%
  full_join(post_eids, by=c("eid"))

dim(all_eids)
which(all_eids[is.na(all_eids)])

all_eids[rowSums(is.na(all_eids)) > 0,]

# read in the outcome data
outcome <- readRDS(outcome_output)
dropped_eids <- c(1077392, 2343986, 2607317, 5422550)

df <- subset(outcome, eid %in% dropped_eids)

head(df)

df %>%
  mutate(epistart = as.Date(epistart, format = "%d/%m/%Y")) %>%
  arrange(eid, epistart) %>%
  group_by(eid) %>%
  mutate(diag = first(epistart)) %>%
  ungroup() %>%
  filter(epistart == diag |is.na(epistart))
  # select(eid, diag_icd10, diag)

head(dropped_outcome)

outcome <- get_first_diag(outcome)
outcome$outcome <- 1
outcome <- unique(outcome[, c("eid", "outcome", "diag")])
length(unique(outcome$eid))
dim(outcome)

## Removing Cancer Cases
icd10_cancer=c(c(paste0("C0", 0:9), paste0("C", 10:14)), paste0("C", 15:26),
               paste0("C", 30:39), paste0("C", 40:41), paste0("C", 43:44),
               paste0("C", 45:49), paste0("C", 50), paste0("C", 51:58),
               paste0("C", 60:63), paste0("C", 64:68), paste0("C", 69:72),
               paste0("C", 73:75), paste0("C", 76:80), paste0("C", 81:96),
               paste0("C", 97), paste0("D0", 0:9), paste0("D", 37:48))
cancer_cases <- get_hes_cases(hes_diag_fname, withdrawn_fname, icd10_cancer)
dim(cancer_cases)
