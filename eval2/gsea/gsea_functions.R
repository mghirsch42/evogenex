library(tidyr)
library(dplyr)
library(clusterProfiler)
library(stringr)

prep_list <- function(data) {
    uniprot_ids <- bitr(data$gene_name, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Mm.eg.db")
    data <- merge(uniprot_ids, data, by.x = "SYMBOL", by.y = "gene_name")
    genelist <- data[,"logFC"]
    names(genelist) <- data[,"UNIPROT"]
    genelist <- sort(genelist, decreasing = TRUE)
}

gse_go <- function(data, save_file) {
    genelist <- prep_list(data)
    try (
        ego <- gseGO(
            geneList = genelist,
            ont = "ALL",
            OrgDb = "org.Mm.eg.db",
            keyType = "UNIPROT",
            pvalueCutoff = 0.05,
            pAdjustMethod = "fdr"
        )
    )
    if (exists("ego")) {
        write.csv(ego, save_file, row.names = FALSE)
    }
    
}


gse_kegg <- function(data, save_file) {
    genelist <- prep_list(data)
    ekegg <- gseKEGG(
        geneList = genelist, 
        organism="mmu", 
        keyType="uniprot",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr"
    )
    write.csv(ekegg, save_file, row.names = FALSE)
}