library(tidyr)
library(dplyr)
library(clusterProfiler)
library(stringr)

gse_kegg <- function(data, save_file) {
    genelist <- data[,"lgfold"]
    names(genelist) <- data[,"UNIPROT"]
    genelist <- sort(genelist, decreasing = TRUE)
    ekegg <- gseKEGG(genelist, organism="mmu", keyType="uniprot")
    write.csv(ekegg, save_file, row.names = FALSE)
}

data_file <- "results/mouse_treatment/ttest_allavg_details.csv"
save_path <- "results/mouse_treatment/gsea_kegg/"

data <- read.csv(data_file)
data["lgfold"] <- log2(data["r.mean.expr"]/data["n.mean.expr"])
uniprot_ids <- bitr(data$gene_name, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Mm.eg.db")
data <- merge(uniprot_ids, data, by.x = "SYMBOL", by.y = "gene_name")

# print(head(genelist))
# results <- gseKEGG(genelist, organism="mmu", keyType="uniprot")
# print(results)
# quit()
for (regime in c("agg", "clade", "X1_4_22", "X3_10_14")) {
# for (regime in c("clade", "X1_4_22", "X3_10_14")) {
# for (regime in c("agg")) {
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
    gse_kegg(query_data, save_file)

    # 2) R > RN; negatively adaptive
    save_file <- paste(save_path, regime, "_R_N.csv", sep = "")
    print(save_file)
    query_data <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] > data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    gse_kegg(query_data, save_file)

    # 3) R < NR; positively adaptive
    save_file <- paste(save_path, regime, "_NR_P.csv", sep = "")
    print(save_file)
    query_data <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] > data[theta_base_col], ]
    gse_kegg(query_data, save_file)

    # 2) R < RN; negatively adaptive
    save_file <- paste(save_path, regime, "_NR_N.csv", sep = "")
    print(save_file)
    query_data <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    gse_kegg(query_data, save_file)

}
