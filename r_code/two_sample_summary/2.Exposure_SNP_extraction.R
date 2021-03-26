# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Install special packages
# devtools::install_github("mrcieu/ieugwasr")

# Read in libraries
library(data.table)
library(ieugwasr)

# Read in the config
config <- read_yaml('../configs/main.yml')

# Read in raw GWAS
SNPs = fread(config$exposure_SNP_fname, stringsAsFactors=FALSE)

#Convert p-values to numbers 
SNPs$pval = as.numeric(SNPs$pval)
plot(SNPs$pval) #checking "uniformity?" 

#Only keep significant SNPs 
sig_SNPs = SNPs[SNPs$pval<5*10^-8,]
dim(sig_SNPs)

#Keep independent SNPs 
ind_SNPs = ieugwasr::ld_clump(sig_SNPs, clump_r2 = 0.001)
dim(ind_SNPs) 

#Save statistically significant independent SNPs 
saveRDS(ind_SNPs, config$exposure_SNP_output)
