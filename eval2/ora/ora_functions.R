library(tidyr)
library(dplyr)
library(clusterProfiler)

kegg_enrich <- function(query_genes, ref_genes, save_file) {
    query_uniprot <- bitr(query_genes$gene_name, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Mm.eg.db")
    ref_uniprot <- bitr(ref_genes$gene_name, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Mm.eg.db")
    print(head(query_uniprot$UNIPROT))
    # quit()
    ego_kegg <- enrichKEGG(
                        gene = query_uniprot$UNIPROT,
                        universe = ref_uniprot$UNIPROT,
                        keyType = "uniprot",
                        organism = "mmu",
                        pvalueCutoff = .05,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = .05
    )
    print(ego_kegg)
    write.csv(ego_kegg, save_file, row.names = FALSE)
}

go_enrich <- function(query_genes, ref_genes, save_file) {
    # print(ref_genes)
    go_all <- enrichGO(
                    gene = query_genes$gene_name,
                    universe = ref_genes$gene_name,
                    OrgDb = "org.Mm.eg.db",
                    keyType = "SYMBOL",
                    ont = "ALL",
                    pvalueCutoff = .05,
                    pAdjustMethod = "fdr",
                    qvalueCutoff = .05,
                    readable = TRUE
                )
    print(go_all)
    write.csv(go_all, save_file, row.names = FALSE)
}