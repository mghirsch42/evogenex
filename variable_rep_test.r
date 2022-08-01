library(tidyverse)
library(ape)
library(ouch)
library(reshape2)
library(nloptr)

tree_file <- "tree_files/resolved/sc-bwes-cons-resolved-10.tree"
regime_file <- "regime_files/resolved/single_resolved.csv"
data_file <- "data/test/tpm_allrep_gene_0-10.csv"
# 
# # Load & modify data as needed
data <- read.csv(data_file)
data <- data[which(data$gene == "ENSMUSG00000000001.4_Gnai3"),]
data <- data[which(data$species != "C2"),]
data$replicate <- paste("R", data$replicate, sep="")


tree_file <- "C:/Users/hirsc/AppData/Local/R/win-library/4.2/EvoGeneX/extdata/drosophila9.newick"
regime_file <- "C:/Users/hirsc/AppData/Local/R/win-library/4.2/EvoGeneX/extdata/regime_global.csv"
data_file <- "C:/Users/hirsc/AppData/Local/R/win-library/4.2/EvoGeneX/extdata/HD_M_FBgn0000008.csv"

# Load & modify data as needed
data <- read.csv(data_file)
rep_data <- data %>% mutate(R5 = R4)
data <- data %>% gather("replicate", "exprval", -species)
rep_data <- rep_data %>% gather("replicate", "exprval", -species)
subset_rep_data <- rep_data
subset_rep_data[sample(which(rep_data$replicate == "R5"), size=3),]$exprval <- NA
subset_rep_data <- subset_rep_data[complete.cases(subset_rep_data),]


run_evogenex <- function(tree_file, regime_file, data) {
  evog <- EvoGeneX()
  evog$setTree(tree_file)
  evog$setRegimes(regime_file)
  results <- evog$fitSlow(data, 0.01, 0.01)
  return(results)
}

run_mcmc_evogenex <- function(tree_file, regime_file, data) {
  evog <- EvoGeneX()
  evog$setTree(tree_file)
  evog$setRegimes(regime_file)
  results <- evog$fitMCMC(data, 0.01, 0.01)
  return(results)
}

# First, run with standard EvoGeneX
library(EvoGeneX)
stand_results <- run_evogenex(tree_file, regime_file, data)
stand_rep_results <- run_evogenex(tree_file, regime_file, rep_data)
old_mcmc_results <- run_mcmc_evogenex(tree_file, regime_file, data)


# Then, run with local EvoGeneX
detach(package:EvoGeneX, unload=TRUE)
source("../EvoGeneX-Lib/Rpackage/R/utils.R")
source("../EvoGeneX-Lib/Rpackage/R/BaseModel.R")
source("../EvoGeneX-Lib/Rpackage/R/EvoGeneX.R")
new_results <- run_evogenex(tree_file, regime_file, data)
new_rep_results <- run_evogenex(tree_file, regime_file, rep_data)
new_subset_rep_results <- run_evogenex(tree_file, regime_file, subset_rep_data)

new_mcmc_results <- run_mcmc_evogenex(tree_file, regime_file, data)

# Compare results
all.equal(stand_results, new_results) # Results with the original data are the same for both versions
all.equal(stand_rep_results, new_rep_results) # Results with replicated data are the same for both versions
all.equal(old_mcmc_results, new_mcmc_results)


# Results with different data will differ from each other
print(new_results)
print(new_rep_results)
print(new_subset_rep_results)
