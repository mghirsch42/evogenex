library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval2/ora/ora_functions.R")

method <- "go" # go or kegg
group <- "negative" # positive, negative, or all
regime <- "agg"
save_file <- paste(
    "results/tpm_allrep2/orig/gene_lists/adpt_", regime,
    "/enrichment/ora_diff_expr_", method,
    "_", group, ".csv",
    sep = ""
)
results_file <- paste(
    "results/tpm_allrep2/orig/gene_lists/adpt_", regime,
    "/differential_expression.csv",
    sep = ""
)

results <- read.csv(results_file)

query_genes <- results %>% filter(fdr_p < 0.05)
ref_genes <- results

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