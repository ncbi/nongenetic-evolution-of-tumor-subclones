.libPaths("~/include/R/") # Needed for cluster, comment for local

library(tidyverse)
library(EvoGeneX)
library(knitr)
library(ape)

##########
### Run EvoGeneX to test for adaptive evolution for all genes in given input file.
### Params:
### Arg 1: File of tree in Newick format with branch lengths
### Arg 2: File of the single regime
### Arg 3: File for the two-regime
### Arg 4: Path to the data file
### Arg 5: Path to save the output file
### See EvoGeneX documentation for requirements for these files
##########

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
newick_file <- args[1]
single_regime_file <- args[2]
two_regime_file <- args[3]
data_file <- args[4]
output_file <- args[5]

# Function to process a single gene
process_single_gene <- function(data_tall) {
  evog$setRegimes(single_regime_file)
  ou1_res <- evog$fit(data_tall, format = "tall", alpha = 0.1, gamma = 0.01)
  
  evog$setRegimes(two_regime_file)
  ou2_res <- evog$fit(data_tall, format = "tall", alpha = 0.1, gamma = 0.01)

  brown_res <- brown$fit(data_tall, format = "tall", gamma = 0.01)
  
  # loglikelihood ratio test EvoGeneX VS replicated Brownian motion
  ou2_vs_bm_pvalue <- 1 - pchisq((ou2_res$loglik - brown_res$loglik) * 2, (ou2_dof - brown_dof))
  ou2_vs_ou1_pvalue <- 1 - pchisq((ou2_res$loglik - brown_res$loglik) * 2, (ou2_dof - ou1_dof))
  results <- tibble("ou1_conv" = ou1_res$optim.diagn$convergence, 
                    "ou1_theta" = ou1_res$theta,
                    "ou1_alpha" = ou1_res$alpha,
                    "ou1_sigma_sq" = ou1_res$sigma.sq,
                    "ou1_gamma" = ou1_res$gamma,
                    "ou1_loglik" = ou1_res$loglik,
                    "ou2_conv" = ou2_res$optim.diagn$convergence, 
                    "ou2_theta" = ou2_res$theta,
                    "ou2_alpha" = ou2_res$alpha,
                    "ou2_sigma_sq" = ou2_res$sigma.sq,
                    "ou2_gamma" = ou2_res$gamma,
                    "ou2_loglik" = ou2_res$loglik,
                    "brown_conv" = brown_res$optim.diagn$convergence, 
                    "brown_theta" = brown_res$theta,
                    "brown_alpha" = brown_res$alpha,
                    "brown_sigma_sq" = brown_res$sigma.sq,
                    "brown_gamma" = brown_res$gamma,
                    "brown_loglik" = brown_res$loglik,
                    "ou2_vs_bm_pvalue" = ou2_vs_bm_pvalue, 
                    "ou2_vs_ou1_pvalue" = ou2_vs_ou1_pvalue)
  return(results)
}

# Setup evogenex and brownian motion models
evog <- EvoGeneX()
evog$setTree(newick_file)

brown <- Brown()
brown$setTree(newick_file)

# Degrees of freedom for each model for stat test
ou1_dof <- 4 # alpha, sigma.sq, theta, gamma
ou2_dof <- 5 # alpha, sigma.sq, theta1, theta2, gamma
brown_dof <- 3 # sigma.sq, theta, gamma
fdr_cutoff <- 0.05

# Run Evogenex over all genes
data_tall <- read.csv(data_file)
results <- (
  data_tall
  %>% group_by(gene)
  %>% summarize(process_single_gene(cur_data()), .groups = "keep")
  %>% mutate(qvalue = p.adjust(max(ou2_vs_bm_pvalue, ou2_vs_ou1_pvalue), method = "fdr"))
  %>% mutate(adaptive = ifelse(qvalue < fdr_cutoff, "adaptive", "not-adaptive"))
)
# print(results)
write.csv(results, output_file, row.names=FALSE)
