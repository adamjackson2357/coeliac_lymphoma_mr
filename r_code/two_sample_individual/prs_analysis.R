# Analysis of the prs

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(ieugwasr)
  library(EBPRS)
  library(ROCR)
  library(ggplot2)
  library(ggpubr)
  library(forestmodel)
  library(xtable)
  library(patchwork)
  library(flextable)
  library(caret)
})
source("../utils/extraction.R")
source("../two_sample_individual/prs.R")
source("../two_sample_individual/model_processing.R")
source("../two_sample_individual/prs_visualisation.R")

# read in the config
config <- read_yaml('../../configs/main.yml')
fields_fname <- config$fields_fname
covars_fname <- config$covars_fname
withdrawn_fname <- config$withdrawn_fname
hes_diag_fname = config$hes_diag_fname
withdrawn_fname <- config$withdrawn_fname
exposure_codes <- config$exposure_codes
outcome_codes <- config$outcome_codes
genotype_fname <- config$genotype_fname
case_covars_output <- config$case_covars_output
gwas_pcs_fname <- config$gwas_pcs_fname
gwas_fname <- config$gwas_fname
clump_threshold <- config$clump_threshold
case_covars_formatted <- config$case_covars_formatted

# get the data
# case_covars <- get_features(covars_fname, hes_diag_fname, withdrawn_fname,
#                             fields_fname, exposure_codes, outcome_codes,
#                             genotype_fname, gwas_fname, clump_threshold,
#                             gwas_pcs_fname)
# saveRDS(case_covars, case_covars_output)
case_covars <- readRDS(case_covars_output)

# do some recoding
case_covars <- feature_recoding(case_covars)

# add in participants with special diets
case_covars <- special_diets(case_covars, covars_fname, withdrawn_fname)

saveRDS(case_covars, case_covars_formatted)

### PRS Histograms

case_covars <- readRDS(case_covars_formatted)
case_covars[,"special_diet"] <- ifelse(case_covars$Gluten_free == 1, "Gluten Free", "Not Gluten Free")

coeliac_hist <- create_hist(case_covars, "Coeliac", "Coeliac")
gluten_hist <- create_hist(case_covars, "special_diet","Special Diet")
nhl_hist <- create_hist(case_covars, "NHL", "NHL")

# Print both side by side
histograms <- ggarrange(
  coeliac_hist, gluten_hist, nhl_hist,ncol=1, nrow=3,
  common.legend = FALSE)

### PRS ROC Curves

## Stage 1: Regress prs against exposure
case_covars <- readRDS(case_covars_output)
case_covars <- gluten_free(case_covars, covars_fname, withdrawn_fname)

# Coeliac logistic regression
exposure_reg <- glm(exposure~prs, data = case_covars, family = binomial(link = "logit"))
exposure_roc <- create_roc_curve(case_covars, exposure_reg, "exposure", "Coeliac")
exposure_prec_rec <- create_prec_rec_curve(case_covars, exposure_reg, "exposure", "Coeliac")

# Gluten free logistic regression
gluten_reg <- glm(Gluten_free~prs, data = case_covars, family = binomial(link = "logit"))
gluten_roc <- create_roc_curve(case_covars, gluten_reg, "Gluten_free", "Gluten Free")
gluten_prec_rec <- create_prec_rec_curve(case_covars, gluten_reg, "Gluten_free", "Gluten Free")

# NHL logistic regression
outcome_reg <- glm(outcome~prs, data = case_covars, family = binomial(link = "logit"))
outcome_roc <- create_roc_curve(case_covars, outcome_reg, "outcome", "NHL")
outcome_prec_rec <- create_prec_rec_curve(case_covars, outcome_reg, "outcome", "NHL")

rocs <- ggarrange(
  exposure_roc, gluten_roc, outcome_roc,ncol=1, nrow=3,
  common.legend = FALSE)
prec_recs <- ggarrange(
  exposure_prec_rec, gluten_prec_rec, outcome_prec_rec, ncol=1, nrow=3,
  common.legend = FALSE)

## Save all the prs plots

prs_plots <- ggarrange(
  histograms, rocs, prec_recs, ncol=3, nrow=1,
  common.legend = FALSE)
ggsave("../../figures/prs_plots.png", prs_plots, width=13, height=11)


### PRS Table

prs_melt <- as.data.table(case_covars) %>%
  mutate(Cohort = "Cohort") %>%
  melt(id.vars=c("eid"),
       measure.vars = c("Cohort", "Coeliac", "NHL", "Sex",
                        "Gluten_free", "Lactose_free", "Low_Calorie", "Vegetarian",
                        "Vegan", "Other", "No_Special_Diet")) %>%
  left_join(select(case_covars, "eid", "PRS"), by="eid")
head(prs_melt)
prs_table <- prs_melt %>%
  group_by(variable, value) %>%
  summarise(n = n(),
            Mean = mean(PRS),
            SD = sd(PRS),
            Min = min(PRS),
            Max = max(PRS))
head(prs_table, 30)
print(xtable(prs_table, type = "latex", digits=4),
      include.rownames=FALSE,
      file = "../../figures/prs_table.tex")
