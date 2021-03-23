library(data.table)
library(openxlsx)

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

### Loading the data

# The fread() function is very efficient for big datasets
# You can choose to read only a few rows at first using nrows (there is 502,538 rows in total this dataset)
mydata=data.frame(fread("../data/ukb26390.csv", nrows=5))

# The first column "eid" is the unique identifier of a participant
print(colnames(mydata)[1:10])
print(mydata[,1])

# The rest of the columns are characteristics of the participants
# (based on questionnaire data, measurements, linkage to other registries, etc).
# They are given as "Field ID"."ID of data collection"."Array ID"

# 1- The field IDs can be linked to explicit descriptions of the variable
# in the data showcase of the UK Biobank website (https://biobank.ndph.ox.ac.uk/showcase/search.cgi)
# where you can type in the field ID or the name of a variable you would like (e.g. you can type "bmi").

# 2- The ID of data collection is an integer: 0 indicates that the information was collected at baseline,
# 1 is at a second time point, and 2 is a third time point. You will see that additional time points are
# available only for some variables, and of these, only a subset participants have data on all time points.
# We recommend to start with data collected at baseline.

# 3- The array ID is available for variables where multiple answers can be collected. For example, 
# see the variable "Qualification" (Field ID 6138) which is about the degree/qualification of participants.
# When filling the questionnaire, participants had to select which of "College/University degree", 
# "A levels", "GCSE", etc they had. One participant can have multiple degrees. The data collected is stored
# in multiple columns for the same variable. For example, for the data on qualifications collected at baseline
# we have 6138.0.0 containing one of the boxes that was ticked by the participant, and 6138.0.1 containing
# data from a second box ticked by the participant, etc. Participants who ticked only one degree will have 
# missing values for the array IDs > 0.

# You can see that there are multiple columns for the qualification variable in the data:
print(colnames(mydata)[grep("6138",colnames(mydata))])

# You can start exploring the variables available for UK Biobank participants using "Browse" in the
# Data Showcase on the UK Biobank website (https://biobank.ndph.ox.ac.uk/showcase/browse.cgi?id=100471&cd=resources) 
# and identify variables that you think would be of interest for your analyses.

# Note that not all variables are available in the dataset you have been provided with, so you might not be 
# able to look into everything you wanted. You can list the field IDs of variables you would like to look into
# as illustrated in the folder "example_extraction". The script in this folder will extract your variables of interest
# from the big "ukb26390.csv" dataset.


### Data recoding

# Categorical variables are coded with some numbers that correspond to each of the categories.
# For these variables, the code is provided in "Codings_Showcase.csv". 
# For example, qualification (field ID 6138) is using data coding 100305 (this information
# can be found on the website). To recode the categories you can use this code:

coding_id="100305"
mycoding=read.csv("../data/Codings_Showcase.csv")
print(head(mycoding))
mycoding_field=mycoding[which(mycoding[,1]==coding_id),]
mycoding_field=mycoding_field[,-1]
rownames(mycoding_field)=mycoding_field[,1]
print(mycoding_field)

# As it is in raw data:
print(mydata$X6138.0.0)

# Recoded categories:
as.character(mycoding_field[as.character(mydata$X6138.0.0),"Meaning"])


### Biomarker measurements

# Biomarker measurements are stored in a different dataset:
bmk=data.frame(fread("../data/ukb27725.csv", nrows=5))

# You can use the "Biomarker_annotation.xlsx" file to identify the biomarkers based on
# their field IDs:
annot=read.xlsx("../data/Biomarker_annotation.xlsx")


### Last but not least - Mandatory step: checking participant consent

# Over the years, some participants have withdrawn consent and asked that the data about them is 
# not used anymore in research. The eids of participants who have withdrawn have been listed in 
# the file "w19266_20200204.csv". Please remember to remove them from your dataset during the data
# preparation stages. 

withdrawn=as.character(read.csv("../data/w19266_20200204.csv")[,1])
print(withdrawn)
