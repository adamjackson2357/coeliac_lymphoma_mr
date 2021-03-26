# Code to extract case eids for a given ICD10 code

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# import libraries
library(data.table)
library(tidyverse)
library(yaml)

# read in the config
config <- read_yaml('../configs/main.yml')

# Define disease codes
icd10 <- config$icd10

# get the biobank data
bb=data.frame(fread(config$covariates_fname))

# remove participants which have withdrawn
withdrawn=read.csv(config$withdrawn_fname)[,1]
bb <- subset(bb, !(eid %in% withdrawn))


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


# save the cases
saveRDS(case_control, config$case_control_output)



#____________________________________________________________________________


#Sanity check (manual extraction of a single ICD code from each time-point)
cases3 <- filter(bb,
                 X40006.0.0 == "C829" |
                   X40006.1.0 == "C829" |
                   X40006.2.0 == "C829" |
                   X40006.3.0 == "C829" |
                   X40006.4.0 == "C829" |
                   X40006.5.0 == "C829" |
                   X40006.6.0 == "C829" |
                   X40006.7.0 == "C829" |
                   X40006.8.0 == "C829" |
                   X40006.9.0 == "C829" |
                   X40006.10.0 == "C829" |
                   X40006.11.0 == "C829" |
                   X40006.12.0 == "C829" |
                   X40006.13.0 == "C829" |
                   X40006.14.0 == "C829" |
                   X40006.15.0 == "C829" |
                   X40006.16.0 == "C829" |
                   X40006.17.0 == "C829" |
                   X40006.18.0 == "C829" |
                   X40006.19.0 == "C829" |
                   X40006.20.0 == "C829" |
                   X40006.21.0 == "C829" |
                   X40006.22.0 == "C829" |
                   X40006.23.0 == "C829" |
                   X40006.24.0 == "C829" |
                   X40006.25.0 == "C829" |
                   X40006.26.0 == "C829" |
                   X40006.27.0 == "C829" |
                   X40006.28.0 == "C829" |
                   X40006.29.0 == "C829" |
                   X40006.30.0 == "C829" |
                   X40006.31.0 == "C829"
)




