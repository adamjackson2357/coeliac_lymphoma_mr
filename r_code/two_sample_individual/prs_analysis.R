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
source("../2sample_individual/prs.R")
source("../2sample_individual/model_processing.R")

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

## Stage 0.1: PRS analysis

case_covars <- readRDS(case_covars_formatted)

# get the legend labels for the histogram
prs_legend_labels <- rbind(
  paste0("Not Diagnosed (n=", count(subset(case_covars, Coeliac == "Not Diagnosed")), ")"),
  paste0("Diagnosed (n=", count(subset(case_covars, Coeliac == "Diagnosed")), ")"))

# run the histogram
prs_hist <- case_covars %>%
  ggplot(aes(x=PRS, colour=Coeliac)) +
  geom_density() +
  theme_minimal() +
  scale_color_manual(labels = prs_legend_labels, values=c("#F8766D", "#00BFC4")) +
  ylab("Density")
prs_hist
ggsave("../../figures/prs_hist.png", prs_hist, width=10, height=5)

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