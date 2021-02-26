# Descriptive Statistics

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
library(yaml)
library(data.table)
library(tidyverse)
library(tableone)
library(forestmodel)

source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

covars_output <- config$covars_output
exposure_output <- config$exposure_output
outcome_output <- config$outcome_output
case_covars_output <- config$case_covars_output
hes_diag_fname <- config$hes_diag_fname
hes_fname <- config$hes_fname
withdrawn_fname <- config$withdrawn_fname
exposure_name <- config$exposure_name
outcome_name <- config$outcome_name

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

outcome <- readRDS(outcome_output)
outcome <- unique(outcome[, c("eid", "diag_icd10")])
outcome %>%
  group_by(diag_icd10) %>%
  count()

# count the number of icd10 codes per participant
num_codes <- table(table(outcome$eid))
c(num_codes, sum = sum(num_codes))


#### Join the exposures, cases and controls ####

# get the covariates
covars <- readRDS(covars_output)
head(covars)

# get exposure first time of diagnosis
exposure <- readRDS(exposure_output)
exposure <- get_first_diag(exposure)
exposure$exposure <- 1
exposure <- unique(exposure[, c("eid", "exposure", "diag")])

# get outcome first time of diagnosis
outcome <- readRDS(outcome_output)
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
old_names <- c("X21022.0.0", "X52.0.0", "X34.0.0", "X31.0.0")
new_names <- c("age_at_recruitment", "birth_month", "birth_year", "sex")
case_covars <- case_covars %>% rename_at(vars(old_names), ~ new_names)

# set dob of birth as the middle of their month of birth
case_covars$dob <- as.Date(paste0(case_covars$birth_year, "-", case_covars$birth_month, "-", "15"), format="%Y-%m-%d")

# get age today
case_covars$age <- as.numeric(difftime(Sys.Date(), case_covars$dob, unit="days")) / 365.25

# age at diagnosis for the exposure
case_covars$age_diag_exposure <- as.numeric(difftime(case_covars$diag_exposure, case_covars$dob, unit="days")) / 365.25

# age at diagnosis for the outcome
case_covars$age_diag_outcome <- as.numeric(difftime(case_covars$diag_outcome, case_covars$dob, unit="days")) / 365.25

# get participants that are immuno-compromised
immuno <- data.frame(fread(hes_diag_fname, select = c("eid", "diag_icd10")))
immuno <- subset(immuno, startsWith(immuno$diag_icd10, "D8"))
case_covars$immuno <- ifelse(case_covars$eid %in% immuno$eid, 1, 0)

# save as an RDS file so that we don't always have to run the above
saveRDS(case_covars, case_covars_output)

## Table One

# read in the case covariates
case_covars <- readRDS(case_covars_output)

head(case_covars)

# convert to cases and controls
case_covars[,outcome_name] <- factor(ifelse(case_covars$outcome == 0, "control", "case"), levels = c("control", "case"))
case_covars[,exposure_name] <- factor(ifelse(case_covars$exposure == 0, "control", "case"), levels = c("control", "case"))
case_covars$sex <- factor(ifelse(case_covars$sex == 0, "female", "male"), levels = c("female", "male"))
case_covars$immuno <- factor(case_covars$immuno)

# vector of the variable names
var_names <- c(exposure_name, "sex", "age", "age_diag_exposure", "age_diag_outcome", "immuno")

# vector of the categorical variable names
cat_vars <- c(exposure_name, "sex", "immuno")

# use table one function
table1 <- CreateTableOne(vars = var_names,
                            strata = outcome_name,
                            data = case_covars,
                            factorVars = cat_vars)
table1_print <- print(table1, showAllLevels = TRUE,
                      quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
table1_print
write.csv(table1_print, "../figures/table1.csv")

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

# add sex
model2 <- glm(outcome~exposure+sex, data=case_covars, family=binomial(link="logit"))
summary(model2)
exp(model2$coefficients)
anova(model1, model2)

# add age
model3 <- glm(outcome~exposure+sex+age, data=case_covars, family=binomial(link="logit"))
summary(model3)
exp(model3$coefficients)
anova(model2, model3)

# add being immuno-compromised
model4 <- glm(outcome~exposure+sex+age+immuno, data=case_covars, family=binomial(link="logit"))
summary(model4)
anova(model3, model4)
round(cbind(exp(cbind(OR = coef(model4),
                      confint(model4))),
            p_value = summary(model4)$coefficients[,4]), 2)
forest_model(model4)

# add in interaction effects
model5 <- glm(outcome~exposure+sex+age+immuno+exposure:sex, data=case_covars, family=binomial(link="logit"))
summary(model5)

