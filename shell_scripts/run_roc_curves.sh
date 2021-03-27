#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

module load anaconda3/personal

cd Group8/Adam/tds/r_code/two_sample_individual/

Rscript prs_roc.R
