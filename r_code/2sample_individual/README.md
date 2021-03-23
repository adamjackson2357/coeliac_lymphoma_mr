Code to run the two-sample individual Mendelian Randomisation Analysis

Generating the Polygenic Risk Score
- prs.R
    - contains functions to calculate the PRS for each eid for a set of SNPs
- prs_analysis.R
    - cohort level analysis of the PRS
    
    
Two Stage Least Squares
- model_processing.R
    - contains functions to create the model input data
- 2sls.R
    - run the 2sls analysis and create figures/tables
