# Create the polygenic risk scores for each snp-eid combination

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
library(yaml)
library(tidyverse)

# read in the config
config <- read_yaml('../configs/main.yml')
genotype_fname <- config$genotype_fname
ivs_output <- config$ivs_output
prs_output <- config$prs_output

# get the prs
prs <- get_prs(genotype_fname, ivs_output)

# save as an rds
saveRDS(prs, prs_output)

## MISSING DATA

case_covars <- readRDS(config$case_covars_output)
geno2 <- data.frame(geno)
geno2$eid <- rownames(geno2)
str(case_covars)
str(geno2)

dim(missing)

missing <- case_covars %>%
  mutate(eid = as.character(eid)) %>%
  select(eid) %>%
  left_join(geno2, by="eid")

rownames(missing) <- missing$eid
missing <- missing[,-1]

# check missing data by column
missing_cols <- apply(missing, MARGIN = 2, FUN = function(y) sum(length(which(is.na(y)))))
missing_cols
barplot(missing_cols)

# check missing data by column
missing_cols2 <- apply(geno, MARGIN = 2, FUN = function(y) sum(length(which(is.na(y)))))
missing_cols2
barplot(missing_cols2)

# how many rows would we have if we only took complete data
geno_count <- nrow(missing)
complete_count <- nrow(na.omit(missing))
null_count <- geno_count - complete_count
null_count

# check the missing data by row
missing_rows <- apply(missing, MARGIN = 1, FUN = function(y) sum(length(which(is.na(y)))))
missing_rows <- missing_rows[missing_rows != 0]
length(missing_rows)
hist(missing_rows)
table(missing_rows)

# Should we impute the missing data?
# Or take the average?

