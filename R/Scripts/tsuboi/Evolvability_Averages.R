# calculating e and c in orthogonal axis defined by ancestral P
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(evolvability)
library(dplyr)
library(ape)
library(evolqg)

# read all species VCV matrices
#setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
#temp = list.files(pattern = "*.csv")
#vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
#names(vcv)  = gsub(".csv", replacement= "", temp)
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/vcv.RData")

# read phylogeny
#filename <- "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
#tree <- read.nexus(filename)
load("~/Desktop/Primatrees.RData")
species <- names(vcv)
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

# remove vcv which are not in the phylogeny
#vcv <- vcv[!names(vcv) %in% setdiff(names(vcv), tree$tip.label)]
#species <- names(vcv)

# get ancestral VCV eigenvectors
all_cov_matrices <- PhyloW(tree, vcv)
ancestral <- getMRCA(tree, tree$tip.label)
anc_vcv <- all_cov_matrices[[as.character(ancestral)]]

# zero first eigenvalue
eig <- eigen(anc_vcv)
D <- diag(eig$values)
V <- eig$vectors
D2 <- D
D2[1,1] <- max(eig$values) * 1e-8
anc_vcv <- V %*% D2 %*% t(V)

# get pcs
pcs <- eigen(anc_vcv)$vectors

# total number of calculations in for
n <- length(species) * ncol(pcs)

# list for results
results_e <- data.frame(
  Species = character(n),
  Eigenvector = integer(n),
  Evolvability = numeric(n),
  Conditional_Evolvability = numeric(n),
  Average_Evolvability = numeric(n),
  Average_Conditional_Evolvability = numeric(n)
)

k <- 1

for(i in seq_along(species)){
  
  covar <- as.matrix(vcv[[species[i]]])
  
  # remove first eingenvalue or put almost zero on it
  eig <- eigen(covar)
  D <- diag(eig$values)
  V <- eig$vectors
  D2 <- D
  D2[1,1] <- max(eig$values) * 1e-8
  covar <- V %*% D2 %*% t(V)
  
  #
  e_mean <- evolvabilityMeans(covar)[1]
  c_mean<- evolvabilityMeans(covar)[4]
  
  for(j in seq_len(ncol(pcs))){
    evo <- evolvabilityBeta(covar, pcs[, j])
    
    results_e[k, ] <- list(
      species[i],
      j,
      evo$e,
      evo$c,
      e_mean,
      c_mean
    )
    
    k <- k + 1
  }
}

#
saveRDS(results_e, "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Averages.RDS")

