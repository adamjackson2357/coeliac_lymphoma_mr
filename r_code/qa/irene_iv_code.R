# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Install special packages
# devtools::install_github("mrcieu/ieugwasr")

### REQUIRED ###

library(data.table)
library(ieugwasr)
library(EBPRS)

# Read in the config
config <- read_yaml('../configs/main.yml')
gwas_fname <- config$gwas_fname
ivs_output <- config$ivs_output

iv <- get_ivs(gwas_fname, 5*10^-8, 0.001)

# save as an rds
saveRDS(iv.clump, ivs_output)

### REQUIRED ###

# Read in libraries
library(data.table)
library(ieugwasr)
library(EBPRS)

source("extraction.R")

# Read in the config
config <- read_yaml('../configs/main.yml')
gwas_fname <- config$gwas_fname

# Read in the gwas data
gwas = read.csv(gwas_fname, stringsAsFactors=FALSE)

# Check out the data
#dim(gwas)
#head(gwas, n=4)

# Rename the columns of your data = must have rsid and pval
colnames(gwas)[1] = "rsid"

# Examine the pvals
gwas$pval = as.numeric(gwas$pval) 
table(gwas$pval<5*10^-6)
table(gwas$pval<10^-8)
table(gwas$pval<5*10^-8)
hist(gwas$pval, breaks=20)
# OOPS, VERY not uniformally distributed...

# Choose your IVs based on a threshold
iv = gwas[gwas$pval<5*10^-8,]
hist(iv$pval, breaks=20)
dim(iv)


# Clump IVs - remove correlated IVs
iv.clump = ieugwasr::ld_clump(iv, clump_r2 = 0.001)
dim(iv.clump) 


# save as an rds
saveRDS(iv.clump, "../data/gwas_ivs.rds")


############################################
############################################
############################################

#EBPRS package calc:

###
# IF you have test data
train <- fread('../data/Dubois.csv') 

# Rename the columns of your data = must have A1, A2, OR, P, SNP
colnames(train)[c(1,2,3)] = c("SNP", "A1", "P")

# path to the plink bfile without extensions 
test <- read_plink("testpath") 

result <- EBPRS(train=traindat, test=plinkfile, N1=4533, N0=10750)
validate(result$S, truey)

####
# IF no test data
train <- fread('../data/Dubois.csv') 

# Rename the columns of your data = must have A1, A2, OR, P, SNP
colnames(train)[c(1,2,3)] = c("SNP", "A1", "P")

# will only provide estimated effect sizes
result = EBPRS(train=traindat, N1=4533, N0=10750) 

#BUT EVERYTHING SAYS TO CLUMP AND THEN FILTER
