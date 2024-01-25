.libPaths("~/include/R/")
library(OUwie)
library(phylobase)
library(phytools)
library(tidyr)
library(dplyr)

##########
### Simulate multiple data replicates under neutral evolution 
### and with the same parameter values
##########

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
	cat("usage: simulate_data.R <tree-file> <outfile> <regime-file> <theta-root> <theta0> <theta1> <sigmasq> <gamma> <alpha> <num-gene> <num-replicate> <replicate-dropout> <regime0> <regime1>\n")
  cat("tree-file: file containing newick string\n")
  cat("outfile: file name for simulated data\n")
  cat("theta-root: expression value at the root\n")
  cat("sigmasq: the sigma squared value\n")
  cat("r: the r value (dispersion parameter)\n")
  cat("num-gene: the number of genes to simulate\n")
  cat("num-replicate: the number of replicates per gene\n")
  cat("replicate-dropout: the rate at which to remove replicates 
      (0 means all genes will have num-replicate replicates)\n")
	quit(status=1)
}

tree_file = args[1] # File containing newick string
outfile = args[2]   # Output file name
theta_root = as.numeric(args[3])  # State value at the root
sigmasq = as.numeric(args[4])
r = as.numeric(args[5])
ngene = as.integer(args[6])
nrep = as.integer(args[7])
rep_drop = as.numeric(args[8])


print(tree_file)
print(outfile)
print(theta_root)
print(sigmasq)
print(r)
print(ngene)
print(nrep)
print(rep_drop)


tree = read.tree(tree_file)
nterm = length(tree$tip.label)
max.age <- max(phytools::nodeHeights(tree))

N = 2*nterm - 1
regimes = read.csv("regime_files/resolved/single_resolved.csv", stringsAsFactors=FALSE)
tree$node.label = regimes$reg[(nterm+1):N]
regimes = regimes[1:nterm,] # take only the leaves
regimes$species = tree$tip.label
regimes = regimes[, c('species', 'regime')]

repnames = paste0('R', 1:nrep)
dat = data.frame()
for (i in 1:ngene) {
  t = OUwie.sim(tree, regimes, theta0 = theta_root, alpha=c(10e-10, 10e-10), sigma.sq = c(sigmasq, sigmasq), theta=c(theta_root,theta_root), root.age=max.age)
  epsilon = matrix(rnbinom(nterm*nrep, size=r, mu=t$X), nterm, nrep) 
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
dat = gather(dat, replicate, exprval, all_of(repnames))
dat <- na.omit(dat)
write.table(dat, file=paste(outfile,"_raw.csv", sep=""), sep=',', quote=F, row.names=F, col.names=T)
dat <- dat %>% mutate(exprval = log2(1+exprval))
write.table(dat, file=outfile, sep=",", quote=F, row.names=F, col.names=T)
