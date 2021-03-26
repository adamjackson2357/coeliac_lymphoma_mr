# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

#Libraries 
library(data.table)


# Read in the config
config <- read_yaml('../configs/main.yml')

#Load in all SNP data for bb 
genetic_data <- readRDS(config$genetic_data_fname)

#Load in independent and significant exposure SNPs
exposure_SNPs <- readRDS(config$exposure_SNP_output)

#Filter all SNP data, keeping only independent and significant SNPs 
genetic_data <- genetic_data[colnames(genetic_data) %in% exposure_SNPs$rsid]

#Make the row names (which are the eids) into a column
setDT(genetic_data, keep.rownames = TRUE)[]

#Name the new column in genetic_data "eid" so we can merge 
colnames(genetic_data)[1] <- "eid"

#Load biobank dataframe with eids and NHL disease status 
case_control <- readRDS(config$case_control_output)

#Join genetic data with biobank particiapants 
case_control_SNP <- merge(case_control, genetic_data, by = "eid")

#Save a datafram with eid, disease status, and statistically significant, indepenent SNPs 
saveRDS(case_control_SNP, config$case_control_SNP_output)
