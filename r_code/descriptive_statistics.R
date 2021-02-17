# Descriptive Statistics

# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# read in libraries
library(data.table)
library(openxlsx)
library(tidyverse)
library(yaml)
library(tableone)
library(forestmodel)
library(forestplot)

source("extraction.R")

# read in the config
config <- read_yaml('../configs/main.yml')

hes_diag_fname <- config$hes_diag_fname
hes_fname <- config$hes_fname
withdrawn_fname <- config$withdrawn_fname
risk_factors <- config$risk_factors
codes <- config$icd10
covariates_fname <- config$covariates_fname
withdrawn_fname <- config$withdrawn_fname
fields_fname <- config$fields_fname

# get the time of first diagnosis
get_first_diag <- function(df) {
  df <- df %>%
    mutate(epistart = as.Date(epistart, format = "%d/%m/%Y")) %>%
    arrange(eid, epistart) %>%
    group_by(eid) %>%
    mutate(diag = first(epistart)) %>%
    ungroup() %>%
    filter(epistart == diag) %>%
    select(eid, diag_icd10, diag)
  return(df)
}

#### Join the coeliacs, cases and controls ####

# get coeliacs first time of diagnosis
coeliacs <- get_cases(hes_diag_fname, hes_fname, withdrawn_fname, risk_factors)
coeliacs <- get_first_diag(coeliacs)
coeliacs$coeliac <- 1
coeliacs <- unique(coeliacs[, c("eid", "coeliac", "diag")])

# get cases first time of diagnosis
cases <- get_cases(hes_diag_fname, hes_fname, withdrawn_fname, codes)
cases <- get_first_diag(cases)
cases$case <- 1
cases <- unique(cases[, c("eid", "case", "diag")])

# outer join coeliacs and cases
coeliacs_cases <- coeliacs %>%
  full_join(cases, by="eid", suffix=c("_coeliac", "_case"))

# get the covariates
fields <- unname(unlist(read.table(fields_fname, header=FALSE)))
fields
covariates <- get_covariates(covariates_fname, withdrawn_fname, fields)

# outer join the cases/coeliacs to the covariates
# enables us to get all the controls
case_covariates <- coeliacs_cases %>%
  full_join(covariates, by="eid") %>%
  mutate(coeliac = if_else(is.na(coeliac), 0, coeliac)) %>%
  mutate(case = if_else(is.na(case), 0, case))

#### Descriptive Statistics

head(case_covariates)

# rename columns
old_names <- c("X21022.0.0", "X52.0.0", "X34.0.0", "X31.0.0")
new_names <- c("age_at_recruitment", "birth_month", "birth_year", "sex")
df <- case_covariates %>% rename_at(vars(old_names), ~ new_names)

# set dob of birth as the middle of their month of birth
df$dob <- as.Date(paste0(df$birth_year, "-", df$birth_month, "-", "15"), format="%Y-%m-%d")

# age at diagnosis for the cases
df$age_diag_case <- as.numeric(difftime(df$diag_case, df$dob, unit="days")) / 365.25

# age at diagnosis for the coeliac cases
df$age_diag_coeliac <- as.numeric(difftime(df$diag_coeliac, df$dob, unit="days")) / 365.25

# age at diagnosis for the coeliac cases
df$age_today <- as.numeric(difftime(Sys.Date(), df$dob, unit="days")) / 365.25

# get participants that are immuno-compromised
immuno <- data.frame(fread(hes_diag_fname, select = c("eid", "diag_icd10")))
immuno <- subset(immuno, startsWith(immuno$diag_icd10, "D8"))
df$immuno <- ifelse(df$eid %in% immuno$eid, 1, 0)

# convert to cases and controls
df$case <- factor(ifelse(df$case == 0, "control", "case"), levels = c("control", "case"))
df$sex <- factor(df$sex)
df$immuno <- factor(df$immuno)

# vector of the variable names
var_names <- c("coeliac", "sex", "age_today", "age_diag_case", "age_diag_coeliac", "immuno")

# vector of the categorical variable names
cat_vars <- c("coeliac", "sex", "immuno")

# use table one function
table1 <- CreateTableOne(vars = var_names,
                            strata = "case",
                            data = df,
                            factorVars = cat_vars)
table1_print <- print(table1, showAllLevels = TRUE,
                      quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
table1_print

# check coeliac immuno compromised cross-tab
addmargins(table(df$coeliac, df$immuno))
addmargins(table(df$case, df$immuno))


#### Logistic regression of coeliacs on cases ####

# table of coeliacs and cases
addmargins(table(df$coeliac, df$case))

# chi square test of the above table
chisq.test(df$coeliac, df$case)

# null model for cases
model0 <- glm(case~1, data=df, family=binomial(link="logit"))

# logistic regression of coeliac on cases
model1 <- glm(case~coeliac, data=df, family=binomial(link="logit"))
summary(model1)
exp(model1$coefficients)
anova(model0, model1)

# add sex
model2 <- glm(case~coeliac+sex, data=df, family=binomial(link="logit"))
summary(model2)
exp(model2$coefficients)
anova(model1, model2)

# add age
model3 <- glm(case~coeliac+sex+age_today, data=df, family=binomial(link="logit"))
summary(model3)
exp(model3$coefficients)
anova(model2, model3)

# add being immuno-compromised
model4 <- glm(case~coeliac+sex+age_today+immuno, data=df, family=binomial(link="logit"))
summary(model4)
anova(model3, model4)
round(cbind(exp(cbind(OR = coef(model4),
                      confint(model4))),
            p_value = summary(model4)$coefficients[,4]), 2)
forestmodel(model4)

# add in interaction effects
model5 <- glm(case~coeliac+sex+age_today+immuno+coeliac:sex, data=df, family=binomial(link="logit"))
summary(model5)

