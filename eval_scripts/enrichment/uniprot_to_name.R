library(clusterProfiler)
library(stringr)
library(org.Mm.eg.db)

column_name <- "geneID" # geneID for ORA, core_enrichment for GSEA
# result_file <- "results/tpm_allrep2/orig/gene_lists/adpt_agg/enrichment/ora/ora_go.csv"
result_file <- "results/mouse_treatment/gsea_go/agg_NR_N.csv"

data <- read.csv(result_file)
if (nrow(data) == 0) {
    print("No rows -- exiting")
    quit()
}

for (i in 1:nrow(data)){
    uniprot_ids <- str_split(data[i, column_name], "/")
    gene_names <- bitr(uniprot_ids[[1]], fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Mm.eg.db")
    data[i, paste(column_name, "_genename", sep="")] <- paste(unique(gene_names$SYMBOL), collapse="/")
}

write.csv(data, result_file, row.names = FALSE)