# Proposed workflow

field_ids.txt: Find relevant field ids and put in the field_ids.txt
->
covariate_extraction.R: run this from the command line to reduce the size of the covariates dataset
->
disease code(s): define the disease codes of interest
->
cases_extraction.R: extract the cases using the disease code(s)
->
control_extraction.R: need to write some code that can select appropriate controls for the cases
->
Analysis
  -> data_exploration.R
