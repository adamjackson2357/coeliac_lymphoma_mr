# model processing functions

remove_cancer_cases <- function(case_covars, hes_diag_fname, withdrawn_fname){
  ## Removing Cancer Cases from the controls
  icd10_cancer=c(c(paste0("C0", 0:9), paste0("C", 10:14)), paste0("C", 15:26),
                 paste0("C", 30:39), paste0("C", 40:41), paste0("C", 43:44),
                 paste0("C", 45:49), paste0("C", 50), paste0("C", 51:58),
                 paste0("C", 60:63), paste0("C", 64:68), paste0("C", 69:72),
                 paste0("C", 73:75), paste0("C", 76:80), paste0("C", 81:96),
                 paste0("C", 97), paste0("D0", 0:9), paste0("D", 37:48))
  cancer_cases <- get_hes_cases(hes_diag_fname, withdrawn_fname, icd10_cancer)
  case_covars <- case_covars[!(case_covars$eid %in% cancer_cases$eid &
                                 case_covars$outcome == 0),]
  return(case_covars)
}


postprocessing <- function(case_covars) {
  # convert to factors
  case_covars$exposure <- as.factor(case_covars$exposure)
  case_covars$outcome <- as.factor(case_covars$outcome)
  case_covars$sex <- as.factor(case_covars$sex)
  
  # convert to age
  case_covars = case_covars %>% rename(age = birth_year) %>% mutate(age = 2020 - age)
  return(case_covars)
}


feature_recoding <- function(case_covars){
  case_covars$sex <- factor(ifelse(case_covars$sex == 1, "Male", "Female"), levels=c("Female", "Male"))
  case_covars$outcome <- factor(ifelse(case_covars$outcome == 1, "Diagnosed", "Not Diagnosed"), levels=c("Not Diagnosed", "Diagnosed"))
  case_covars$exposure <- factor(ifelse(case_covars$exposure == 1, "Diagnosed", "Not Diagnosed"), levels=c("Not Diagnosed", "Diagnosed"))
  old_names <- c("age", "sex", "outcome", "exposure", "prs")
  new_names <- c("Age", "Sex", "NHL", "Coeliac", "PRS")
  case_covars <- case_covars %>% rename_at(vars(old_names), ~ new_names)
  return(case_covars)
}


gluten_free <- function(case_covars, covars_fname, withdrawn_fname){
  # add in people that have different special diets
  diet <- get_covariates(covars_fname, withdrawn_fname, "20086")
  special_diet <- data.frame(
    eid = diet$eid,
    Gluten_free = ifelse(diet$eid %in% filter_all(diet, any_vars(. %in% 8))$eid, 1, 0))
  case_covars <- case_covars %>%
    left_join(special_diet, by=c("eid"))
  case_covars["Gluten_free"] <- lapply(case_covars["Gluten_free"], as.factor)
  return(case_covars)
}


special_diets <- function(case_covars, covars_fname, withdrawn_fname){
  # add in people that have different special diets
  diet <- get_covariates(covars_fname, withdrawn_fname, "20086")
  special_diet <- data.frame(
    eid = diet$eid,
    Gluten_free = ifelse(diet$eid %in% filter_all(diet, any_vars(. %in% 8))$eid, 1, 0),
    Lactose_free = ifelse(diet$eid %in% filter_all(diet, any_vars(. %in% 9))$eid, 1, 0),
    Low_Calorie = ifelse(diet$eid %in% filter_all(diet, any_vars(. %in% 10))$eid, 1, 0),
    Vegetarian = ifelse(diet$eid %in% filter_all(diet, any_vars(. %in% 11))$eid, 1, 0),
    Vegan = ifelse(diet$eid %in% filter_all(diet, any_vars(. %in% 12))$eid, 1, 0),
    Other = ifelse(diet$eid %in% filter_all(diet, any_vars(. %in% 13))$eid, 1, 0))
  special_diet <- special_diet %>%
    mutate(No_Special_Diet = rowSums(.[2:7]),
           No_Special_Diet = ifelse(No_Special_Diet == 0, 1, 0))
  case_covars <- case_covars %>%
    left_join(special_diet, by=c("eid"))
  
  special_diet_cols <- c("Gluten_free", "Lactose_free", "Low_Calorie", "Vegetarian",
                         "Vegan", "Other", "No_Special_Diet")
  case_covars[special_diet_cols] <- lapply(case_covars[special_diet_cols], as.factor)
  return(case_covars)
}


get_features <- function(covars_fname, hes_diag_fname, withdrawn_fname,
                         fields_fname, exposure_codes, outcome_codes,
                         genotype_fname, gwas_fname, clump_threshold,
                         gwas_pcs_fname){
  case_covars <- get_case_covars(covars_fname, hes_diag_fname, withdrawn_fname,
                                 fields_fname, exposure_codes, outcome_codes)
  
  # recoding
  old_names <- c("X34.0.0", "X31.0.0")
  new_names <- c("birth_year", "sex")
  case_covars <- case_covars %>% rename_at(vars(old_names), ~ new_names)
  
  # get the prs and join to the case covariates dataframe
  prs <- get_prs(genotype_fname, gwas_fname, 5*10^-8, clump_threshold)
  case_covars$prs <- prs[match(case_covars$eid, names(prs))]
  
  # get the first 10 gwas pcs and join
  gwas_pcs <- readRDS(gwas_pcs_fname)[,1:11]
  case_covars <- inner_join(case_covars, gwas_pcs, by=c("eid"))
  table(case_covars$outcome)
  
  # drop any remaining nas
  case_covars <- na.omit(case_covars)
  
  # remove cancer cases from the controls
  case_covars <- remove_cancer_cases(case_covars, hes_diag_fname, withdrawn_fname)
  
  # do some postprocessing
  case_covars <- postprocessing(case_covars)
  
  return(case_covars)
}