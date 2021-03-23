# Check the number of lymphoma cases

# Clear variables and set the path
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(yaml)
library(dplyr)
library(data.table)

source("extraction.R")

config <- read_yaml('../configs/main.yml')
hes_diag_fname <- config$hes_diag_fname
hes_fname <- config$hes_fname
withdrawn_fname <- config$withdrawn_fname

# get the hes_diag data
hes_diag <- data.frame(fread(hes_diag_fname))

# get the hes data
hes <- data.frame(fread(hes_fname))

head(hes_diag)
head(hes)

# inner join
cases <- hes_diag %>%
  inner_join(hes, by = c("eid", "ins_index"))
head(cases)
