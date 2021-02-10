#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

module load anaconda3/personal

cd /rds/general/project/hda_students_data/live/Group8/adam/r_code

Rscript cases_extraction.R
