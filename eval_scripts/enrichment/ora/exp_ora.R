library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval_scripts/enrichment/ora/ora_functions.R")

regime <- "3_10_14"
method <- "kegg" # go or kegg
group <- "negative" # positive, negative, or anything for all
save_file <- paste("results/tpm_allrep2/orig/gene_lists/adpt_",
                    regime, "/enrichment/entrezid/", regime, "_", group,
                    ".csv", sep="")
results_file <- paste("results/tpm_allrep2/orig/gene_lists/adpt_",
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