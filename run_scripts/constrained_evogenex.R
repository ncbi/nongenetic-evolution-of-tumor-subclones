.libPaths("~/include/R/") # Needed for cluster

library(tidyverse)
library(EvoGeneX)
library(knitr)
library(ape)

##########
### Run EvoGeneX to test for constrained evolution for all genes in given input file.
### Params:
### Arg 1: File of tree in Newick format with branch lengths
### Arg 2: File of the single regime
### Arg 3: Path to the data file
### Arg 4: Path to save the output
### See EvoGeneX documentation for requirements for these files
##########

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
newick_file <- args[1]
regime_file <- args[2]
data_file <- args[3]
output_file <- args[4]

# Function to process a single gene
process_single_gene <- function(data_tall) {
  ou_res <- evog$fit(data_tall, format = "tall", alpha = 0.1, gamma = 0.01)
  brown_res <- brown$fit(data_tall, format = "tall", gamma = 0.01)
  # loglikelihood ratio test EvoGeneX VS replicated Brownian motion
  pvalue <- 1 - pchisq((ou_res$loglik - brown_res$loglik) * 2, (ou_dof - brown_dof))
}

# Setup evogenex and brownian motion models
evog <- EvoGeneX()
evog$setTree(newick_file)
evog$setRegimes(regime_file)

brown <- Brown()
brown$setTree(newick_file)

# Setup stat test parameters
ou_dof <- 4 # alpha, sigma.sq, theta, gamma
brown_dof <- 3 # sigma.sq, theta, gamma
fdr_cutoff <- 0.05

# Run Evogenex over all genes
data_tall <- read.csv(data_file)
results <- (
  data_tall
  %>% group_by(gene)
  %>% summarize(pvalue = process_single_gene(cur_data()), .groups = "keep")
  # %>% ungroup()
  %>% mutate(qvalue = p.adjust(pvalue, method = "fdr"))
  %>% mutate(constrained_vs_neutral = ifelse(qvalue < fdr_cutoff, "constrained", "neutral"))
)

write.csv(results, output_file, row.names=FALSE)
