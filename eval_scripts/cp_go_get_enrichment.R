library(tidyr)
library(dplyr)
library(clusterProfiler)

args <- commandArgs(trailingOnly=TRUE)
query_file <- args[1]
ref_file <- args[2]
save_file <- args[3]

query_genes <- read.csv(query_file)
ref_genes <- read.csv(ref_file)

query_genes <- query_genes %>% separate(ensemble_id, c("ensemble_id", NA))
ref_genes <- ref_genes %>% separate(ensemble_id, c("ensemble_id", NA))

ego_bp <- enrichGO(
                    gene = query_genes$ensemble_id,
                    universe = ref_genes$ensemble_id,
                    OrgDb = "org.Mm.eg.db",
                    keyType = "ENSEMBL",
                    ont = "BP",
                    pvalueCutoff = .05,
                    qvalueCutoff = .05,
                    readable = TRUE
                )
ego_cc <- enrichGO(
                    gene = query_genes$ensemble_id,
                    universe = ref_genes$ensemble_id,
                    OrgDb = "org.Mm.eg.db",
                    keyType = "ENSEMBL",
                    ont = "CC",
                    pvalueCutoff = .05,
                    qvalueCutoff = .05,
                    readable = TRUE
                )
                
ego_mf <- enrichGO(
                    gene = query_genes$ensemble_id,
                    universe = ref_genes$ensemble_id,
                    OrgDb = "org.Mm.eg.db",
                    keyType = "ENSEMBL",
                    ont = "MF",
                    pvalueCutoff = .05,
                    qvalueCutoff = .05,
                    readable = TRUE
                )

# png("figures/cp_enrichment/scaled_agg/bp_barplot.png")
# barplot(ego_bp)
# png("figures/cp_enrichment/scaled_agg/cc_barplot.png")
# barplot(ego_cc)
# png("figures/cp_enrichment/scaled_agg/mf_barplot.png")
# barplot(ego_mf)

# dev.off()

ego_bp_df <- as.data.frame(ego_bp)
ego_bp_df$type <- "bio"
ego_cc_df <- as.data.frame(ego_cc)
ego_cc_df$type <- "cell"
ego_mf_df <- as.data.frame(ego_mf)
ego_mf_df$type <- "mol"
ego_all <- rbind(ego_bp_df, ego_cc_df, ego_mf_df)

write.csv(ego_all, save_file, row.names=FALSE)