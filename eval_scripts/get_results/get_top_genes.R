library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
result_file <- args[2]

data <- read.csv(data_file)

top_genes <- bind_rows(
    data 
        %>% filter(logFC >= 0 & qvalue <= 0.01)
        %>% arrange(desc(abs(logFC)))
        %>% head(10),
    data
        %>% filter(logFC <= 0 & qvalue <= 0.01)
        %>% arrange(desc(abs(logFC)))
        %>% head(10)
)

top_genes <- subset(top_genes, select=c("ensemble_id", "gene_name", "ou2_theta", "ou2_theta_base", "theta_diff", "rel_diff", "qvalue", "logFC"))
# write.table(top_genes, result_file, row.names=FALSE)

# print(head(data %>% arrange(desc(data["logFC"]))))
