.libPaths("~/include/R/")
library(tidyverse)
library(dplyr)

##########
### This takes a file from get_gene_results.py and adds columns for qvalue, relative difference of
### theta values, and the log fold change of theta values and save it back to the same file.
### We use R for this to be consistent with the previous EvoGeneX study
##########

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]

data <- read.csv(data_file)
data <- (data 
            %>% group_by(gene_name)
            %>% mutate(qvalue = p.adjust(max(ou2_vs_bm_pvalue, ou2_vs_ou1_pvalue), method = "fdr"))
            %>% mutate(rel_diff = ((ou2_theta - ou2_theta_base) / ou2_theta_base))
            %>% mutate(logFC = log2(ou2_theta/ou2_theta_base))
)

write.csv(data, data_file, row.names=FALSE)

