.libPaths("~/include/R/") # Needed for cluster, comment for local

library(tidyverse)
library(EvoGeneX)
library(knitr)
library(ape)

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
newick_file <- args[1]
single_regime_file <- args[2]
three_regime_file <- args[3]
data_file <- args[4]
output_file <- args[5]

# Function to process a single gene
process_single_gene <- function(data_tall) {
  evog$setRegimes(single_regime_file)
  ou1_res <- evog$fit(data_tall, format = "tall", alpha = 0.1, gamma = 0.01)
  
  evog$setRegimes(three_regime_file)
  ou3_res <- evog$fit(data_tall, format = "tall", alpha = 0.1, gamma = 0.01)

  brown_res <- brown$fit(data_tall, format = "tall", gamma = 0.01)
  
  # loglikelihood ratio test EvoGeneX VS replicated Brownian motion
  ou3_vs_bm_pvalue <- 1 - pchisq((ou3_res$loglik - brown_res$loglik) * 2, (ou3_dof - brown_dof))
  ou3_vs_ou1_pvalue <- 1 - pchisq((ou3_res$loglik - brown_res$loglik) * 2, (ou3_dof - ou1_dof))
  results <- tibble("ou1_conv" = ou1_res$optim.diagn$convergence, 
                    "ou1_theta" = ou1_res$theta,
                    "ou1_alpha" = ou1_res$alpha,
                    "ou1_sigma_sq" = ou1_res$sigma.sq,
                    "ou1_gamma" = ou1_res$gamma,
                    "ou1_loglik" = ou1_res$loglik,
                    "ou3_conv" = ou3_res$optim.diagn$convergence, 
                    "ou3_theta_agg_nr" = ou3_res$theta["aggressive_nonresponder"],
                    "ou3_theta_agg_r" = ou3_res$theta["aggressive_responder"],
                    "ou3_theta_other" = ou3_res$theta["other"],
                    "ou3_alpha" = ou3_res$alpha,
                    "ou3_sigma_sq" = ou3_res$sigma.sq,
                    "ou3_gamma" = ou3_res$gamma,
                    "ou3_loglik" = ou3_res$loglik,
                    "brown_conv" = brown_res$optim.diagn$convergence, 
                    "brown_theta" = brown_res$theta,
                    "brown_alpha" = brown_res$alpha,
                    "brown_sigma_sq" = brown_res$sigma.sq,
                    "brown_gamma" = brown_res$gamma,
                    "brown_loglik" = brown_res$loglik,
                    "ou3_vs_bm_pvalue" = ou3_vs_bm_pvalue, 
                    "ou3_vs_ou1_pvalue" = ou3_vs_ou1_pvalue)
  return(results)
}

# Setup evogenex and brownian motion models
evog <- EvoGeneX()
evog$setTree(newick_file)

brown <- Brown()
brown$setTree(newick_file)

# Setup stat test parameters
ou1_dof <- 4 # alpha, sigma.sq, theta, gamma
ou3_dof <- 6 # alpha, sigma.sq, theta1, theta2, theta3, gamma
brown_dof <- 3 # sigma.sq, theta, gamma
fdr_cutoff <- 0.05

# Run Evogenex over all genes
data_tall <- read.csv(data_file)
results <- (
  data_tall
  %>% group_by(gene)
  %>% summarize(process_single_gene(cur_data()), .groups = "keep")
  %>% mutate(qvalue = p.adjust(max(ou3_vs_bm_pvalue, ou3_vs_ou1_pvalue), method = "fdr"))
  %>% mutate(adaptive = ifelse(qvalue < fdr_cutoff, "adaptive", "not-adaptive"))
)
# print(results)
write.csv(results, output_file, row.names=FALSE)
