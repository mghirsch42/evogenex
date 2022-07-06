# .libPaths("~/include/R/") # Needed for cluster, comment for local

library(tidyverse)
library(EvoGeneX)
library(knitr)
library(ape)

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
newick_file <- args[1]
single_regime_file <- args[2]
two_regime_file <- args[3]
data_file <- args[4]
output_file <- args[5]

# Function to process a single gene
process_single_gene <- function(data_tall) {
  evog$setRegimes(single_regime_file)
  ou1_res <- evog$fit(data_tall, format = "tall", alpha = 0.1, gamma = 0.01)
  
  evog$setRegimes(two_regime_file)
  ou2_res <- evog$fit(data_tall, format = "tall", alpha = 0.1, gamma = 0.01)

  brown_res <- brown$fit(data_tall, format = "tall", gamma = 0.01)
  
  # loglikelihood ratio test EvoGeneX VS replicated Brownian motion
  ou2_vs_bm_pvalue <- 1 - pchisq((ou2_res$loglik - brown_res$loglik) * 2, (ou2_dof - brown_dof))
  ou2_vs_ou1_pvalue <- 1 - pchisq((ou2_res$loglik - brown_res$loglik) * 2, (ou2_dof - ou1_dof))
  results <- tibble("ou2_vs_bm_pvalue" = ou2_vs_bm_pvalue, "ou2_vs_ou1_pvalue" = ou2_vs_ou1_pvalue)
  return(results)
}

# Setup evogenex and brownian motion models
evog <- EvoGeneX()
evog$setTree(newick_file)

brown <- Brown()
brown$setTree(newick_file)

# Setup stat test parameters
ou1_dof <- 4 # alpha, sigma.sq, theta, gamma
ou2_dof <- 5 # alpha, sigma.sq, theta1, theta2, gamma
brown_dof <- 3 # sigma.sq, theta, gamma
fdr_cutoff <- 0.05

# Run Evogenex over all genes
data_tall <- read.csv(data_file)
results <- (
  data_tall
  %>% group_by(gene)
  %>% summarize(process_single_gene(cur_data()), .groups = "keep")
  %>% mutate(adaptive = ifelse(max(ou2_vs_bm_pvalue, ou2_vs_ou1_pvalue) < fdr_cutoff, "adaptive", "not-adaptive"))
)
# print(results)
write.csv(results, output_file, row.names=FALSE)
