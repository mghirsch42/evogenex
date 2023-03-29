library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval_scripts/enrichment/ora/ora_functions.R")

method <- "kegg" # go or kegg
# Define what group to look at based on their theta values. 
# 2 means that regime has the max value, 
# 1 means that regime has the middle value, 
# 0 means that regime has the min value.
# Exactly one regime should be assigned to each value.
agg_nr <- 2 # 0, 1, 2
agg_r <- 0 # 0, 1, 2
other <- 1 # 0, 1, 2

save_file <- paste(#"results/tpm_allrep2/orig/gene_lists/adpt_three/enrichment/ora_",
                   # method, 
                    "_agg_nr_", agg_nr,
                    "_agg_r_", agg_r,
                    "_other_", other,
                    ".csv",
                    sep = "")
results_file <- "results/tpm_allrep2/orig/gene_lists/adpt_three/all_results.csv"
results <- read.csv(results_file)

# Filter based on theta conditions
if (agg_nr > agg_r) {
    results <- results %>% filter(ou3_theta_agg_nr > ou3_theta_agg_r)
} else {
    results <- results %>% filter(ou3_theta_agg_nr < ou3_theta_agg_r)
}
if (agg_nr > other) {
    results <- results %>% filter(ou3_theta_agg_nr > ou3_theta_other)
} else {
    results <- results %>% filter(ou3_theta_agg_nr < ou3_theta_other)
}
if (agg_r > other) {
    results <- results %>% filter(ou3_theta_agg_r > ou3_theta_other)
} else {
    results <- results %>% filter(ou3_theta_agg_r < ou3_theta_other)
}

# Query genes are those that are adaptive;
query_genes <- (
    results 
    %>% filter(adaptive == "adaptive")
    %>% separate(ensemble_id, c("ensemble_id", NA))
)
ref_genes <- (
    results
    %>% separate(ensemble_id, c("ensemble_id", NA))
)

print("Number of query genes:")
print(nrow(query_genes))

# write.csv(query_genes$gene_name, file=save_file, row.names=FALSE, quote=FALSE)
# quit()
if (method == "go") {
    go_enrich(query_genes, ref_genes, save_file)
} else {
    kegg_enrich(query_genes, ref_genes, save_file)
}