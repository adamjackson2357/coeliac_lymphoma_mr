# TDS Project

## Introduction

Group 8 analysis for TDS Module.

Using Mendelian Randomisation to assess the causative effects of coeliac disease on
non-hodgkins lymphoma.

## Quick start

1. Clone the project from [https://github.com/adamjackson2357/coeliac_lymphoma_mr]
2. Create the mr_env environment
```
module load anaconda3/personal
conda env create configs/mr_env.yml
source activate mr_env
```

## Project Structure

```bash
├── configs
│   ├── field_ids.txt
│   ├── main.yml
│   └── README.md
├── data
├── documentation
│   ├── README_original.txt
│   └── workflow.md
├── r_code
│   ├── cases_extraction.R
│   ├── covariate_extraction.R
│   ├── data_exploration.R
│   ├── example_extraction
│   │   ├── List_field_ids_to_extract.txt
│   │   └── script_covariates_extraction.R
│   ├── function_example.R
│   ├── getting_started.R
│   ├── README.md
│   └── recoding_disease.R
└── README.md
```

To create this tree, navigate in the terminal to the project root and type `tree`
Copy and paste the output into within the pairs of triple backticks.
This will change quite a lot and is worth updating as we add new files
Files in the /data directory shouldn't be uploaded to git

## Prerequisites

r_packages:
- devtools
- yaml
- data.table
- tidyverse
- tableone -
- forestmodel
- ieugwasr
- EBPRS
- ROCR
- ivprobit

UK biobank data

## Documentation

The documentation for this project is this README.

See `/documentation/workflow.md` for the proposed workflow of different scripts to run

See `/documentation/README_original.txt` for the recommended steps on getting started.

## Support team

- Ahmed Abdulaal
- Aimee Tham
- Irene Sugurova
- Adam Jackson
