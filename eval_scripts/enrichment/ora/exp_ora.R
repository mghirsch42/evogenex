library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval2/ora/ora_functions.R")

method <- "kegg" # go or kegg
group <- "positive" # positive, negative, or all
save_file <- "results/tpm_allrep2/orig/gene_lists/adpt_agg/enrichment/ora_kegg_positive.csv"
results_file <- "results/tpm_allrep2/orig/gene_lists/adpt_agg/all_results.csv"
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