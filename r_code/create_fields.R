# Data exploration code

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# read in libraries
library(tidyverse)
library(yaml)

# read in the config
config <- read_yaml('../configs/main.yml')

filter_field <- "selected_20210127"

# read in the data dictionary
data_dictionary <- read_csv(config$data_dict_fname)

# filter for the select column
fields <- subset(data_dictionary, !is.na(data_dictionary[,filter_field]))$field_id

# save the field_id 
write_tsv(fields, paste("../configs/field_ids_", sub(".*_", "", filter_field), ".txt", sep=""))
