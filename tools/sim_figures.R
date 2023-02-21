#!/usr/bin/env Rscript
library(OUwie)
library(phylobase)
library(phytools)

# Simulate data with adaptive evolution with the same alpha and sigma squared

tree_file = "tree_files/resolved/sc-bwes-cons-resolved-10.tree" # File containing newick string
outfile = "fake"   # Output file name
regime_file = "regime_files/resolved/agg_resolved.csv" # Regime file
theta_root = 150  # State value at the root
theta_ratio = 1.1
sigmasq = 1.5
gamma = 1
alpha = 1
ngene = 100
nrep = 8
regime0 = "nonaggressive"
regime1 = "aggressive"

tree = read.tree(tree_file)
nterm = length(tree$tip.label)
max.age <- max(phytools::nodeHeights(tree))

N = 2*nterm - 1
regimes = read.csv(regime_file, stringsAsFactors = FALSE)
regimes[regimes == regime0] = 1 # relabel the regimes so OUSim will apply the correct theta
regimes[regimes == regime1] = 2
tree$node.label = regimes$reg[(nterm+1):N]  # label the tree with the regimes

regimes = regimes[1:nterm,] # take only the leaves
regimes$species = tree$tip.label
regimes = regimes[, c('species', 'regime')]

repnames = paste0('R', 1:nrep)

theta0=(2*theta_root)/(1+theta_ratio)
theta1=(2*theta_root*theta_ratio) / (1+theta_ratio)
print(theta0)
print(theta1)

dat = data.frame()

for (i in 1:ngene) {
  t = OUwie.sim(tree, regimes, theta0 = theta_root, alpha = c(alpha,alpha), 
                sigma.sq = c(sigmasq, sigmasq), theta = c(theta0, theta1), 
                root.age = max.age, get.all = TRUE)
  print(t)
  quit()
  epsilon = matrix(rnbinom(nterm*nrep, size=100, mu=t$X), nterm, nrep) 
  for (s in 1:nterm) {
    for (j in 1:nrep) {
      if (runif(1, 0, 1) < rep_drop) {
        epsilon[s, j] = NA
      }
    }
  }
  Y = data.frame(epsilon)
  names(Y) = repnames
  Y$species = t$Genus_species
  gene = i
  dat = rbind(dat, data.frame(gene, Y))
}
write.table(dat[c('gene', 'species', repnames)], file=outfile, sep='\t', quote=F, row.names=F, col.names=T)
