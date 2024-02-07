### EPIMERGER
#
# Created by: Roc Farriol
#    XX XX 2022
# 
# Modified by: Victor Montal
#    7 Feb 2024

# GOAL: MERGE EPITOPES from different methods into a single dataframe

# --
## Libraries
# --
rm(list = ls())
library(tidyr, quietly = T, warn.conflicts = F)
library(tibble, quietly = T, warn.conflicts = F)
library(dplyr, quietly = T, warn.conflicts = F)
library(data.table, quietly = T, warn.conflicts = F)
library(argparser, quietly = T, warn.conflicts = F)


# --
## SCRIPT PARSER
# --
# Create a parser
p <- arg_parser("EPIMERGER")
p <- add_argument(p, "--path", help= "Path to the brewpitope project", type="character")
argv <- parse_args(p)


# --
## Define defaults 
# --
# Lineal
ibepi_lin <- paste0(argv$path,"/C_epixtractor/bepipred_results_extracted.csv")
iabc_lin <- paste0(argv$path,"/C_epixtractor/abcpred_results_extracted.csv")
iepivec_lin <- paste0(argv$path,"/C_epixtractor/epitopevec_results_extracted.csv")

# Structural
iseppa_conf <- paste0(argv$path,"/C_epixtractor/seppa_results_extracted.csv")
idisco_conf <- paste0(argv$path,"/C_epixtractor/discotope_results_extracted.csv") 
ibepi_conf <- paste0(argv$path,"/C_epixtractor/bepitope_conf_results_extracted.csv")
iserendip_conf <- paste0(argv$path,"/C_epixtractor/serendipce_results_extracted.csv")

# --
## Load different dataframes
# --
epi_lin_df <- read.csv(file = ibepi_lin , sep = ",")
epivec_lin_df <- read.csv(file = iepivec_lin , sep = ",")
seppa_conf_df <- read.csv(file = iseppa_conf , sep = ",")
disco_conf_df <- read.csv(file = idisco_conf , sep = ",")
bepi_conf_df <- read.csv(file = ibepi_conf , sep = ",")
serendip_conf_df <- read.csv(file = iserendip_conf , sep = ",")
  
# --
## Merge tables
# --
merged <- bind_rows(epi_lin_df, epivec_lin_df,seppa_conf_df,disco_conf_df,bepi_conf_df,serendip_conf_df)

# --
## Export
# --
ofile <- paste0(argv$path,"/D_epimerger/merged.csv")
write.table(merged, file = ofile, row.names = F, quote = F, sep = ";")
print(paste("Find your output merged file at: ",ofile))
