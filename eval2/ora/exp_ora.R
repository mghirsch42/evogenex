library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval2/ora/ora_functions.R")

method <- "kegg" # go or kegg
save_file <- "results/tpm_allrep2/orig/gene_lists/adpt_3_10_14/ora_kegg.csv"
results_file <- "results/tpm_allrep2/orig/gene_lists/adpt_3_10_14/all_results.csv"
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

if (method == "go") {
    go_enrich(query_genes, ref_genes, save_file)
} else {
    kegg_enrich(query_genes, ref_genes, save_file)
}