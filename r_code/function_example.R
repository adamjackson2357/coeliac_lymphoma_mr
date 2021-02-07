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
