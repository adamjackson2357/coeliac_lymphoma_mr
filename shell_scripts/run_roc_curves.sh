#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=20gb

module load anaconda3/personal
source activate mr_env

cd /rds/general/user/aj1520/home/Group8/Adam/tds/r_code/two_sample_individual/

Rscript prs_roc.R
