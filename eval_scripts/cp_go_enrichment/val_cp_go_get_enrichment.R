library(tidyr)
library(dplyr)
library(clusterProfiler)

bio_enrich <- function(query_genes, ref_genes) {
    ego_bp <- enrichGO(
                        gene = query_genes,
                        universe = ref_genes,
                        OrgDb = "org.Mm.eg.db",
                        keyType = "SYMBOL",
                        ont = "BP",
                        pvalueCutoff = .05,
                        qvalueCutoff = .05,
                        readable = TRUE
                    )
}

cell_enrich <- function(query_genes, ref_genes) {
    ego_cc <- enrichGO(
                        gene = query_genes,
                        universe = ref_genes,
                        OrgDb = "org.Mm.eg.db",
                        keyType = "SYMBOL",
                        ont = "CC",
                        pvalueCutoff = .05,
                        qvalueCutoff = .05,
                        readable = TRUE
                    )
}

mol_enrich <- function(query_genes, ref_genes) {
    ego_mf <- enrichGO(
                        gene = query_genes,
                        universe = ref_genes,
                        OrgDb = "org.Mm.eg.db",
                        keyType = "SYMBOL",
                        ont = "MF",
                        pvalueCutoff = .05,
                        qvalueCutoff = .05,
                        readable = TRUE
                    )
}

all_enrich <- function(query_genes, ref_genes, save_file) {
    go_bio <- bio_enrich(query_genes$gene_name)
    go_cell <- cell_enrich(query_genes$gene_name)
    go_mol <- mol_enrich(query_genes$gene_name)

    go_all <- data.frame()
    bio_df <- as.data.frame(go_bio)
    if (nrow(bio_df) > 0) {
        bio_df$type <- "bio"
        go_all <- rbind(go_all, bio_df)
    }
    cell_df <- as.data.frame(go_cell)
    if (nrow(cell_df) > 0) {
        cell_df$type <- "cell"
        go_all <- rbind(go_all, cell_df)
    }
    mol_df <- as.data.frame(go_mol)
    if (nrow(mol_df) > 0) {
        mol_df$type <- "mol"
        go_all <- rbind(go_all, mol_df)
    }
    write.csv(go_all, save_file, row.names = FALSE)
}

args <- commandArgs(trailingOnly=TRUE)
# query_file <- args[1]
data_file <- "results/mouse_treatment/ttest_allavg_details.csv"
# save_file <- args[3]
save_path <- "results/mouse_treatment/enrichment/"

data <- read.csv(data_file)
ref_genes <- data$gene_name
# for (regime in c("agg", "clade", "X1_4_22", "X3_10_14")) {
for (regime in c("clade", "X1_4_22", "X3_10_14")) {
# for (regime in c("agg")) {
    print(regime)
    query_col <- paste(regime, "q.value", sep = ".")
    theta_col <- paste(regime, "ou2.theta", sep = ".")
    theta_base_col <- paste(regime, "ou2.theta.base", sep = ".")

    # 1) R > NR; positively adaptive
    save_file <- paste(save_path, regime, "_cp_go_R_P.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] > data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] > data[theta_base_col], ]
    print(nrow(query_genes))
    all_enrich(query_genes, ref_genes, save_file)

    # 2) R > RN; negatively adaptive
    save_file <- paste(save_path, regime, "_cp_go_R_N.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] > data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    print(nrow(query_genes))
    all_enrich(query_genes, ref_genes, save_file)


    # 3) R < NR; positively adaptive
    save_file <- paste(save_path, regime, "_cp_go_NR_P.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] > data[theta_base_col], ]
    print(nrow(query_genes))
    all_enrich(query_genes, ref_genes, save_file)

    # 2) R < RN; negatively adaptive
    save_file <- paste(save_path, regime, "_cp_go_NR_N.csv", sep = "")
    print(save_file)
    query_genes <- data[data["fdr.p.value"] < 0.05
                                & data["r.mean.expr"] < data["n.mean.expr"]
                                & data[query_col] < 0.05
                                & data[theta_col] < data[theta_base_col], ]
    print(nrow(query_genes))
    all_enrich(query_genes, ref_genes, save_file)
}
