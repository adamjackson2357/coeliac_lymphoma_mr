# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(xtable)
library(ieugwasr)
library(dplyr)
source("../2sample_individual/prs.R")

variants <- read.csv("../../data/Marq.csv")
head(variants)
ivs <- get_ivs("../../data/Marq.csv", 5*10^-8, 0.001)
head(ivs)

ivs$EA <- ivs$MA

old_names <- c("rsid", "EA","MA", "MAF", "pval", "OR")
new_names <- c("rsid", "Effect Allele", "Minor Allele", "MAF", "p-value", "Odds Ratio")
ivs <- ivs %>% rename_at(vars(old_names), ~ new_names)

ivs <- ivs[,new_names]

head(ivs)

print(xtable(ivs, type = "latex", digits = 3, display=c("s","s","s", "s","s","g","s")),
      include.rownames=FALSE,
      math.style.exponents = TRUE,
      file = "../../figures/instrument_variables.tex")
