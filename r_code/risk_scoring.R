# Create the risk scores for each snp-eid combination

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in libraries
library(yaml)

# read in the config
config <- read_yaml('../configs/main.yml')
genotype_fname <- config$genotype_fname
rs_output <- config$rs_output

## Read in the data

# read in the snp data
geno <- readRDS(genotype_fname)

# select snps
snps <- c("rs7705526",  "rs6920364",  "rs11780471", "rs4236709")
or <- c(0.5, 1.5, 1, 2)
risk_allele <- c("minor", "minor", "major", "major")

# subset the genotype data to only the snps of interest
geno <- geno[, snps]

# convert fom raw to numeric
for (i in 1:length(snps)) {geno[, snps[i]] <- as.numeric(geno[, snps[i]])}

# convert to a matrix
geno <- as.matrix(geno)

## Genotype Re-coding

# convert all 0s to NAs
geno[geno == 0] <- NA

# if it's a minor allele then convert 3 to 0
# if it's a major allele then convert 3 to 1, and 1 to 0
pre_process <- function(x, risk_allele) {
  if (risk_allele == "minor"){
    x[x == 3] <- 0
  } else {
    x[x == 1] <- 0
    x[x == 3] <- 1
  }
  return(x)
}

# apply to each column
for (i in 1: length(snps)) {geno[, i] <- pre_process(geno[,i], risk_allele[i])}

## Imputation

# Should we impute the missing data

## Create the risk score

# multiply each SNP by the odds ratio to get the risk for that snp
rs <- t(t(geno) * or)

# save as an rds
saveRDS(rs, rs_output)
