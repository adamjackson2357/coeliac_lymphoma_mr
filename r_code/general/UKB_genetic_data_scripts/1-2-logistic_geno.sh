#!/bin/bash
#PBS -l walltime=30:00:00
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -N gwas_geno_m1
#PBS -q pqeph
#PBS -J 1-22

cd /work/bbodinie/UK_Biobank/Covid19_severity/Scripts
module load plink

geno_path=/rds/general/project/uk-biobank-2018/live/reference/sdata_12032018
fam_path=/work/bbodinie/UK_Biobank/Covid19_severity/Data
data_path=/work/bbodinie/UK_Biobank/Covid19_severity/Data
results_path=/work/bbodinie/UK_Biobank/Covid19_severity/Results

plink --bfile $geno_path/ukb_gen_chr$PBS_ARRAY_INDEX \
--fam $fam_path/ukb_geno_m1.fam \
--covar $data_path/confounders.txt keep-pheno-on-missing-cov --no-const-covar \
--maf 0.01 \
--ci 0.95 \
--covar-name BatchBatch_b002-BatchUKBiLEVEAX_b9,PC1-PC10 \
--hide-covar \
--logistic --out $results_path/logistic_geno_m1_chr$PBS_ARRAY_INDEX
