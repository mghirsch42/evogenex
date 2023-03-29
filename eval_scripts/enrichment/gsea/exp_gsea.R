library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval_scripts/enrichment/gsea/gsea_functions.R")

method <- "kegg" # go or kegg
group <- "positive" # positive, negative, or all
save_file <- "results/tpm_allrep2/orig/gene_lists/adpt_1_4_22/enrichment/gsea_kegg_positive.csv"
results_file <- "results/tpm_allrep2/orig/gene_lists/adpt_1_4_22/all_results.csv"
results <- read.csv(results_file)

data <- results %>% filter(adaptive == "adaptive")

if (group == "positive") {
    data <- data %>% filter(logFC > 0)
} else if (group == "negative") {
    data <- data %>% filter(logFC < 0)
}

if (method == "kegg") {
    gse_kegg(data, save_file)
} else {
    gse_go(data, save_file)
}
