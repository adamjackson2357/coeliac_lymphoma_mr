# Code to extract case eids for a given ICD10 code

# import libraries
library(data.table)
library(tidyverse)
library(yaml)

# read in the config
config <- read_yaml('../configs/main.yml')

# Define disease codes
codes <- config$codes

# get the hospital statistics data
hes=data.frame(fread(config$hes_fname))

# remove participants which have withdrawn
withdrawn=read.csv(config$withdrawn_fname)[,1]
hes <- subset(hes, !(eid %in% withdrawn))

# get unique eids back for a given disease code
cases <- hes %>%
  filter(diag_icd10 %in% codes) %>%
  distinct(eid)

# save the cases
saveRDS(cases, config$cases_output)
