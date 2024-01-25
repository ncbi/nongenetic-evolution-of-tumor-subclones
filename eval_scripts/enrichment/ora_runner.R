library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval_scripts/enrichment/ora_functions.R")

##########
### Script to run over-representation enrichment analysis using clusterProfiler
### on the results from EvoGeneX.
### Options to run with KEGG or GO, although only KEGG was used in the results of this study.
##########

regime <- "har"
method <- "kegg" # go or kegg
group <- "negative" # positive or negative
save_file <- paste("results/tpm/gene_lists/adpt_",
                    regime, "/enrichment/", regime, "_", group,
                    ".csv", sep="")
results_file <- paste("results/tpm/gene_lists/adpt_",
                       regime, "/all_results.csv", sep="")
results <- read.csv(results_file)

query_genes <- (
    results 
    %>% filter(adaptive == "adaptive")
    %>% separate(ensemble_id, c("ensemble_id", NA))
)
ref_genes <- (
    results
    %>% separate(ensemble_id, c("ensemble_id", NA))
)

if (group == "positive") {
    query_genes <- query_genes %>% filter(logFC > 0)
} else if (group == "negative") {
    query_genes <- query_genes %>% filter(logFC < 0)
}

if (method == "go") {
    go_enrich(query_genes, ref_genes, save_file)
} else {
    kegg_enrich(query_genes, ref_genes, save_file)
}