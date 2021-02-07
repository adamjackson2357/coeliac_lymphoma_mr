# TDS Project

## Introduction

UK Biobank analysis for TDS Module.

## Quick start

1. Clone the project from [git_lab_link]

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
- tidyverse
- data.table
- yaml
- openxlsx

UK biobank data

## Useful Info

README.md
- READMEs in the .md (markdown) format will render really nicely on git
- You can find them in each directory
- Please update them (including this one!) as we go!
- I think we have a lecture on markdown later on in TDS

.gitignore 
- Stops git from pushing some files to the remote (bad practice to push sensitive data)
- It's 'hidden' in the terminal, use `ls -a` to see that it's there

Documentation
- Found this standardised documentation for R functions called Roxygen2
- https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html

Coding Tips!
- Hopefully these don't come across as patronising! Just a few things I've picked up over the last few years.
- Try to avoid capital letters as much as possible. E.g. in variable names, filenames etc. This will make everything much easier in the long run. This README is an unfortunate contradiction.
- Use '_' instead of spaces in strings.
- A consistent date format is really useful; I'd suggest `YYYYMMDD` as it makes sorting easy. E.g. "20210813", 13th August 2021

## Documentation

The documentation for this project is this README.

See `/documentation/workflow.md` for the proposed workflow of different scripts to run

See `/documentation/README_original.txt` for the recommended steps on getting started.

## Support team

-- Ahmed Abdulaal
-- Aimee Tham
-- Irene Sugurova
-- Adam Jackson
