library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval_scripts/enrichment/ora/ora_functions.R")

method = "kegg"
data_file <- "results/mouse_treatment/ttest_allavg_details.csv"
save_path <- "results/mouse_treatment/ora_kegg_entrez/"

data <- read.csv(data_file)
ref_genes <- data
for (regime in c("agg", "clade", "X1_4_22", "X3_10_14")) {
# for (regime in c("clade", "X1_4_22", "X3_10_14")) {
# for (regime in c("agg")) {
    print(regime)
    query_col <- paste(regime, "q.value", sep = ".")
    theta_col <- paste(regime, "ou2.theta", sep = ".")
    theta_base_col <- paste(regime, "ou2.theta.base", sep = ".")

    # 1) R > NR; positively adaptive
    save_file <- paste(save_path, regime, "_R_P.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] > data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] > data[theta_base_col], ]
    print(nrow(query_genes))
    if (method == "kegg") {
        kegg_enrich(query_genes, ref_genes, save_file)
    }
    else {
        go_enrich(query_genes, ref_genes, save_file)
    }
    # 2) R > RN; negatively adaptive
    save_file <- paste(save_path, regime, "_R_N.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] > data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    print(nrow(query_genes))
    if (method == "kegg") {
        kegg_enrich(query_genes, ref_genes, save_file)
    }
    else {
        go_enrich(query_genes, ref_genes, save_file)
    }

    # 3) R < NR; positively adaptive
    save_file <- paste(save_path, regime, "_NR_P.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] > data[theta_base_col], ]
    print(nrow(query_genes))
    if (method == "kegg") {
        kegg_enrich(query_genes, ref_genes, save_file)
    }
    else {
        go_enrich(query_genes, ref_genes, save_file)
    }
    # 2) R < RN; negatively adaptive
    save_file <- paste(save_path, regime, "_NR_N.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    print(nrow(query_genes))
    if (method == "kegg") {
        kegg_enrich(query_genes, ref_genes, save_file)
    }
    else {
        go_enrich(query_genes, ref_genes, save_file)
    }
}
