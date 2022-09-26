library(biomaRt)
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
print(nrow(query_genes))
print(nrow(ref_genes))
# Convert ensemble genes to uniprot
print("converting ids")
ensemble <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
query_entrez <- getBM(attributes = c("ensembl_gene_id", "uniprot_gn_id"),
                      filters = "ensembl_gene_id",
                      values = query_genes$ensemble_id,
                      mart = ensemble)
ref_entrez <- getBM(attributes = "uniprot_gn_id",
                    filters = "ensembl_gene_id",
                    values = ref_genes$ensemble_id,
                    mart = ensemble)

# Remove nulls - genes that couldn't be mapped
query_entrez <- query_entrez %>% drop_na()
ref_entrez <- ref_entrez %>% drop_na()
print(nrow(query_entrez))
print(nrow(ref_entrez))

# Get enrichment
print("getting enrichment")
ekegg <- enrichKEGG(
                    gene = query_entrez$uniprot_gn_id,
                    universe = ref_entrez$uniprot_gn_id,
                    organism = "mmu",
                    keyType = "uniprot",
                    pvalueCutoff = .05,
                    qvalueCutoff = .05
                )
print(ekegg)
ekegg_df <- as.data.frame(ekegg)
write.csv(ekegg, save_file)