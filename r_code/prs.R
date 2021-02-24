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

## Read in the data

# read in the data
geno <- readRDS(genotype_fname)
ivs <- readRDS(ivs_output)

# re-order the genotype data
geno <- geno[, ivs$rsid]

# convert to a numeric matrix while preserving the eids
eids <- rownames(geno)
geno <- as.matrix(sapply(geno, as.numeric))
rownames(geno) <- eids

# convert all 0s to NAs
geno[geno == 0] <- NA

# convert to minor homozygote to 2
# the heterozygote to 1
# and the major homozygote to 0
recode_alleles <- function(x) {abs(x - 3)}

# demonstrate
x <- c(1, 2, 3)
recode_alleles(x)

# apply to the geno dataset
geno <- apply(geno, MARGIN = 2, FUN = recode_alleles)

# multiply each SNP by the log odds ratio
rs <- t(t(geno) * log(ivs$OR))
head(rs)

# additive sum and exponentiate to get the prs
prs <- exp(rowSums(rs, na.rm = TRUE))

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

