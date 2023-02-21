library(tidyr)
library(dplyr)
library(clusterProfiler)
source("eval2/gsea/gsea_functions.R")

# print(sessionInfo())
# quit()
data_file <- "results/mouse_treatment/ttest_allavg_details.csv"
save_path <- "results/mouse_treatment/gsea_go/"
method <- "go"

data <- read.csv(data_file)
data["logFC"] <- log2(data["r.mean.expr"]/data["n.mean.expr"])

for (regime in c("agg", "clade", "X1_4_22", "X3_10_14")) {
    print(regime)
    query_col <- paste(regime, "q.value", sep = ".")
    theta_col <- paste(regime, "ou2.theta", sep = ".")
    theta_base_col <- paste(regime, "ou2.theta.base", sep = ".")
    
    # # 1) R > NR; positively adaptive
    save_file <- paste(save_path, regime, "_R_P.csv", sep = "")
    print(save_file)
    query_data <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] > data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] > data[theta_base_col], ]
    if (method == "kegg") {
        gse_kegg(query_data, save_file)
    }
    else {
        print(nrow(query_data))
        gse_go(query_data, save_file)
    }

    # 2) R > RN; negatively adaptive
    save_file <- paste(save_path, regime, "_R_N.csv", sep = "")
    print(save_file)
    query_data <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] > data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    if (method == "kegg") {
        gse_kegg(query_data, save_file)
    }
    else {
        print(nrow(query_data))
        gse_go(query_data, save_file)
    }
    

    # 3) R < NR; positively adaptive
    save_file <- paste(save_path, regime, "_NR_P.csv", sep = "")
    print(save_file)
    query_data <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] > data[theta_base_col], ]
    if (method == "kegg") {
        gse_kegg(query_data, save_file)
    }
    else {
        print(nrow(query_data))
        gse_go(query_data, save_file)
    }

    # 2) R < RN; negatively adaptive
    save_file <- paste(save_path, regime, "_NR_N.csv", sep = "")
    print(save_file)
    query_data <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    if (method == "kegg") {
        gse_kegg(query_data, save_file)
    }
    else {
        print(nrow(query_data))
        gse_go(query_data, save_file)
    }
}
