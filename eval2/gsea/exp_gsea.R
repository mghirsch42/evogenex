library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval2/gsea/gsea_functions.R")

method <- "kegg" # go or kegg
save_file <- "results/tpm_allrep2/orig/gene_lists/adpt_clade/gsea_kegg.csv"
results_file <- "results/tpm_allrep2/orig/gene_lists/adpt_clade/all_results.csv"
results <- read.csv(results_file)

data <- results %>% filter(adaptive == "adaptive")

if (method == "kegg") {
    gse_kegg(data, save_file)
} else {
    gse_go(data, save_file)

}
