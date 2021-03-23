# Descriptive Statistics

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(tidyverse)
  library(tableone)
  library(xtable)
  library(forestmodel)
})
source("../utils/extraction.R")

# read in the config
config <- read_yaml('../../configs/main.yml')

covars_output <- config$covars_output
exposure_output <- config$exposure_output
outcome_output <- config$outcome_output
case_covars_output <- config$case_covars_output
hes_diag_fname <- config$hes_diag_fname
hes_fname <- config$hes_fname
withdrawn_fname <- config$withdrawn_fname
exposure_name <- config$exposure_name
outcome_name <- config$outcome_name
gwas_pcs_fname <- config$gwas_pcs_fname
outcome_codes <- config$outcome_codes
exposure_codes <- config$exposure_codes
fields <- config$fields_fname
covars_fname <- config$covars_fname

column_ids <- data.frame(fread(covars_fname, nrows=1))

head(column_ids)

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

#### Get Statistics for each ICD10 code ####

outcome <- get_cases_duration(hes_diag_fname, hes_fname, withdrawn_fname, outcome_codes)
outcome <- unique(outcome[, c("eid", "diag_icd10")])
outcome %>%
  group_by(diag_icd10) %>%
  count()

# count the number of icd10 codes per participant
num_codes <- table(table(outcome$eid))
c(num_codes, sum = sum(num_codes))


#### Join the exposures, cases and controls ####

fields <- unname(unlist(read.table(config$fields_fname, header=FALSE)))

# get the covariates
covars <- get_covariates(covars_fname, withdrawn_fname, fields)
head(covars)

# get exposure first time of diagnosis
exposure <- get_cases_duration(hes_diag_fname, hes_fname, withdrawn_fname, exposure_codes)
exposure <- get_first_diag(exposure)
exposure$exposure <- 1
exposure <- unique(exposure[, c("eid", "exposure", "diag")])

# get outcome first time of diagnosis
outcome <- get_cases_duration(hes_diag_fname, hes_fname, withdrawn_fname, outcome_codes)
outcome <- get_first_diag(outcome)
outcome$outcome <- 1
outcome <- unique(outcome[, c("eid", "outcome", "diag")])

# outer join the exposure/outcome to the covariates
# enables us to get all the controls
case_covars <- covars %>%
  left_join(exposure, by="eid") %>%
  left_join(outcome, by="eid", suffix=c("_exposure", "_outcome")) %>%
  mutate(exposure = if_else(is.na(exposure), 0, exposure)) %>%
  mutate(outcome = if_else(is.na(outcome), 0, outcome))

#### Descriptive Statistics Pre-processing

# rename columns
old_names <- c("X34.0.0", "X31.0.0")
new_names <- c("birth_year", "sex")
case_covars <- case_covars %>% rename_at(vars(old_names), ~ new_names)

# set dob of birth as the middle of their month of birth
case_covars$dob <- as.Date(paste0(case_covars$birth_year, "-", case_covars$birth_month, "-", "15"), format="%Y-%m-%d")

# get age today
Sys.Date()
case_covars$age <- as.numeric(difftime(Sys.Date(), case_covars$dob, unit="days"))

# age at diagnosis for the exposure
case_covars$age_diag_exposure <- as.numeric(difftime(case_covars$diag_exposure, case_covars$dob, unit="days")) / 365.25

# age at diagnosis for the outcome
case_covars$age_diag_outcome <- as.numeric(difftime(case_covars$diag_outcome, case_covars$dob, unit="days")) / 365.25

# save as an RDS file so that we don't always have to run the above
saveRDS(case_covars, case_covars_output)

## Additional Data

# read in the case covariates
case_covars <- readRDS(case_covars_output)

# get participants that are immuno-compromised
icd10 <- data.frame(fread(hes_diag_fname, select = c("eid", "diag_icd10")))
immuno <- subset(icd10, startsWith(icd10$diag_icd10, "D8"))
case_covars$immuno <- ifelse(case_covars$eid %in% immuno$eid, 1, 0)

# get hiv cases
hiv <- rbind(subset(icd10, startsWith(icd10$diag_icd10, "B20")),
             subset(icd10, startsWith(icd10$diag_icd10, "B21")),
             subset(icd10, startsWith(icd10$diag_icd10, "B22")),
             subset(icd10, startsWith(icd10$diag_icd10, "B23")),
             subset(icd10, startsWith(icd10$diag_icd10, "B24")))
hiv <- unique(hiv$eid)
case_covars$hiv <- ifelse(case_covars$eid %in% hiv, 1, 0)

## Table One
head(case_covars)

# convert to cases and controls
case_covars[,outcome_name] <- factor(ifelse(case_covars$outcome == 0, "control", "case"), levels = c("control", "case"))
case_covars[,exposure_name] <- factor(ifelse(case_covars$exposure == 0, "control", "case"), levels = c("control", "case"))
case_covars$sex <- factor(ifelse(case_covars$sex == 0, "female", "male"), levels = c("female", "male"))
case_covars$immuno <- factor(case_covars$immuno)

# vector of the variable names
var_names <- c(exposure_name, "sex", "age", "immuno", "hiv")

# vector of the categorical variable names
cat_vars <- c(exposure_name, "sex", "immuno", "hiv")

# use table one function
table1 <- CreateTableOne(vars = var_names,
                            strata = outcome_name,
                            data = case_covars,
                            factorVars = cat_vars)
table1_print <- print(table1, showAllLevels = FALSE,
                      quote = FALSE, noSpaces = TRUE, test = FALSE)

table1_print
print(xtable(table1_print, type = "latex"), file = "../../figures/table_1.tex")

head(case_covars)
case_covars %>% 
  group_by(outcome) %>%
  summarise(min=min(age),
            max=max(age),
            median=median(age),
            IQR = IQR(age))


## Venn Diagram

# subset participants that are cases, coeliacs or immunocompromised
overlap <- subset(case_covars, outcome == "lymphoma" | exposure == "coeliac" | immuno == 1)
overlap <- overlap[,c("outcome", "exposure", "immuno")]

overlap_agg <- overlap %>%
  group_by(outcome, exposure, immuno) %>%
  count() %>%
  arrange(outcome, exposure, immuno)

head(overlap_agg, 10)

#### Logistic regression of exposures on outcomes ####

case_covars <- readRDS(case_covars_output)
case_covars <- head(case_covars, 1000)
table(case_covars$outcome)

case_covars$sex <- factor(ifelse(case_covars$sex == 1, "Male", "Female"), levels=c("Female", "Male"))
case_covars$outcome <- factor(ifelse(case_covars$outcome == 1, "Diagnosed", "Not Diagnosed"), levels=c("Not Diagnosed", "Diagnosed"))
case_covars$exposure <- factor(ifelse(case_covars$exposure == 1, "Diagnosed", "Not Diagnosed"), levels=c("Not Diagnosed", "Diagnosed"))
case_covars$immuno <- factor(ifelse(case_covars$immuno == 1, "Diagnosed", "Not Diagnosed"), levels=c("Not Diagnosed", "Diagnosed"))
case_covars$hiv <- factor(ifelse(case_covars$hiv == 1, "Diagnosed", "Not Diagnosed"), levels=c("Not Diagnosed", "Diagnosed"))
  
old_names <- c("age", "sex", "outcome", "exposure", "immuno", "hiv")
new_names <- c("Age", "Sex", "NHL", "Coeliac", "Immunologically_Deficient", "HIV")
case_covars <- case_covars %>% rename_at(vars(old_names), ~ new_names)

# get the different function strings
pc_cols <- names(readRDS(gwas_pcs_fname))[2:11]
outcome_str <- "NHL ~ "
covars_str <-  paste(c("Sex", "Age", "Coeliac", "Immunologically_Deficient", "HIV", pc_cols), collapse='+')

mini_case_covars <- head(case_covars, 50000)
model <- glm(as.formula(c(outcome_str, covars_str)), data=case_covars, family=binomial(link="logit"))
summary(model)
covars_list = c("Sex", "Age", "Coeliac", "Immunologically_Deficient", "HIV")
format_list <- forest_model_format_options(banded = FALSE)
obs_log_reg <- forest_model(model, factor_separate_line=TRUE,
                            covariates=covars_list, format_options = format_list)
ggsave(plot = obs_log_reg, filename="../../figures/obs_forestplot.png", width=10, height=5)

# table of exposures and outcomes
addmargins(table(case_covars$exposure, case_covars$outcome))

# chi square test of the above table
chisq.test(case_covars$exposure, case_covars$outcome)

# null model for outcomes
model0 <- glm(outcome~1, data=case_covars, family=binomial(link="logit"))

# logistic regression of exposure on outcomes
model1 <- glm(outcome~exposure, data=case_covars, family=binomial(link="logit"))
summary(model1)
exp(model1$coefficients)
anova(model0, model1)

# get the different function strings
pc_cols <- names(readRDS(gwas_pcs_fname))[2:11]
outcome_str <- "outcome ~ "
covars_str <-  paste(c("age", "sex", pc_cols), collapse='+')
covars_exp_str <- paste(c(covars_str, "exposure"), collapse='+')

# add in immuno and hiv
covars_immuno_str <-  paste(c("age", "sex", "immuno", "hiv", pc_cols), collapse='+')
covars_immuno_exp_str <- paste(c(covars_immuno_str, "exposure"), collapse='+')

model2 <- glm(as.formula(c(outcome_str, covars_str)), data=case_covars, family=binomial(link="logit"))
summary(model2)
exp(model2$coefficients)

model3 <- glm(as.formula(c(outcome_str, covars_exp_str)), data=case_covars, family=binomial(link="logit"))
summary(model3)
exp(model3$coefficients)
anova(model2, model3)
forest_model(model3, reference)
?forest_model

model4 <- glm(as.formula(c(outcome_str, covars_immuno_str)), data=case_covars, family=binomial(link="logit"))
summary(model4)
exp(model4$coefficients)

model5 <- glm(as.formula(c(outcome_str, covars_immuno_exp_str)), data=case_covars, family=binomial(link="logit"))
summary(model5)
exp(model5$coefficients)
anova(model4, model5)

## Forest Plot

model5

forest_model(model5)

