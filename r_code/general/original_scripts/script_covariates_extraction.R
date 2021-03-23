# As the entire UK Biobank dataset is quite big, it can be useful to create a copy with only the covariates of interest
# (i.e. discarding any irrelevant variable for your analyses)
# The code below is creating a copy containing only the variables of interest as defined in
# a list of field IDs ("List_field_ids_to_extract.txt")
# The field IDs of variables can be found on the UK Biobank website (in the data showcase)
# Note: You need to change the path where the data is copied

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(data.table)

# Loading the data
mydata=data.frame(fread("../../data/ukb26390.csv", nrows=5))
myfields=unname(unlist(read.table("List_field_ids_to_extract.txt", header=FALSE)))

# Extracting the column ids 
column_id=grep("eid", colnames(mydata))
found_fieldids=NULL
for (k in 1:length(myfields)){
mygrep=grep(paste0("X",myfields[k],"."), fixed=TRUE, colnames(mydata))
if (length(mygrep)>0){
found_fieldids=c(found_fieldids, myfields[k])
}
column_id=c(column_id, mygrep)
}

# Extracting required columns from dataset
extracted=data.frame(fread("../../data/ukb26390.csv", select=column_id))
saveRDS(extracted, "covid_individual_covariates.rds")

