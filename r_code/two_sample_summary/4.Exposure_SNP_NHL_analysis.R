# Clear variables and set the path
dev.off()
rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

# Read in the config
config <- read_yaml('../configs/main.yml')

#Load case/control + exposure SNPs data
case_control_SNP <- readRDS(config$case_control_SNP_output)

#Convert the SNP data to numeric data-type
case_control_SNP[,16:38] <- sapply(case_control_SNP[,16:38],as.numeric)

#Convert "0" genotypes to NA 
case_control_SNP[,16:38][case_control_SNP[,16:38] == 0] <- NA

#Recoding alleles so that 2 = double minor allele, 1 = heterozygous minor, and 0 = double major
recode_SNPs <- function(x) {
  ifelse(x == 3, 0, 
         ifelse(x == 2, 1, 2))
}

case_control_SNP[,16:38] <- apply(case_control_SNP[,16:38], MARGIN = 2, function(x) recode_SNPs(x))

#Save recoded allele data 
saveRDS(case_control_SNP, config$case_control_SNP_recoded_output)

colnames(case_control_SNP)

#Initialise some matrices for collecting required data 
SNP <- NULL
beta <- NULL
se <- NULL
pval <- NULL

#Point to the principle component columns of the dataframe
PC <- paste(colnames(case_control_SNP[6:15]), collapse='+')

#Look through the SNPs and perform univariable logistic regression
for (i in 16:length(colnames(case_control_SNP))) {
  
  #Save SNP name  
  SNP <- rbind(SNP, colnames(case_control_SNP[i]))
  
  formula <- "outcome ~"
  formula <- c(formula, paste(colnames(case_control_SNP[i]), "+age+sex+"))
  formula <- c(formula, PC)
  formula <- as.formula(formula)
  
  #Run univariable logistic regression model
  model <- glm(formula, family=binomial (link=logit),
               data = case_control_SNP)
  
  #Save beta values 
  beta <- rbind(beta, model$coefficients[[2]])
  
  #Save standard errors 
  se <- rbind(se, sqrt(diag(vcov(model)))[[2]])
  
  #Save p-values 
  pval <- rbind(pval, coef(summary(model))[,4][[2]])
}

#Two Sample MR
library(devtools)
library(TwoSampleMR)
library(MRInstruments)

#Load preparatory outcomes data (contains SNP, effect allele, other allele, and allele frequency) 
outcomes_data <- read.csv("../data/outcomes_data.csv")

#Create data-frame with SNP, betas, ses and p-values from previous analysis 
gwas <-as.data.frame(cbind(SNP, beta, se, pval))

#Give columns default names expected by TwoSampleMR package
x <- c("SNP", "beta", "se", "pval")
colnames(gwas) <-x

#Sanity check data frames
gwas
outcomes_data

#Merge into one - reading for formatting
gwas <- merge(gwas, outcomes_data, by="SNP")

#Defining a factor -> numeric function (not actually available in R)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#Coercing data into correct types 
gwas$beta <- as.numeric.factor(gwas$beta)
gwas$se <- as.numeric.factor(gwas$se)
gwas$pval <- as.numeric.factor(gwas$pval)
gwas$effect_allele <- as.character(gwas$effect_allele)
gwas$other_allele <- as.character(gwas$other_allele)
gwas

#Format the outcome data for MR analysis 
outcome_data <- format_data(gwas, type="outcome")
outcome_data

#Exposure data 
exposure_data <- read.csv("../data/exposure_data.csv")
exposure_data <- format_data(exposure_data, type="exposure")
exposure_data

#harmonize (just in case)
dat <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)

#MR analysis 
res <- mr(dat)
res

#Hetero
mr_heterogeneity(dat)

#Pleiotropy 
mr_pleiotropy_test(dat)

#Funnel plot 
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Leave one out + plot
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

#Single SNP analysis 
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

#Outlier detection - MR PRESSO 
run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)


qnorm(6.69E-09)
