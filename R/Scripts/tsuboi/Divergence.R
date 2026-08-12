# Calculate divergence between Ps
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ape)
library(evolqg)

# function to calculate diff
responseDiff <- function(P1, P2, beta){
  
  D <- P1 - P2
  
  sqrt(as.numeric(t(beta) %*% D %*% D %*% beta))
}

# scale trace of P matrix
scale.trace <- function(P, target.trace){
  P * target.trace / sum(diag(P))
}

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

# calculating the divergence between species
Pavg <- anc_vcv

# generate 1000 vectors
ntraits <- nrow(Pavg)
nrep <- 1000
beta <- matrix(rnorm(ntraits * nrep), ntraits, nrep)
beta <- apply(beta, 2, function(x) x / sqrt(sum(x^2)))

target.trace <- sum(diag(as.matrix(Pavg)))

results_d <- data.frame()

for(sp in names(vcv)){
  
  P <- as.matrix(vcv[[sp]])
  Pavg <- as.matrix(Pavg)
  
  ## sem padronizar o trace
  d.raw <- apply(beta, 2, responseDiff,
                 P1 = P,
                 P2 = Pavg)
  
  ## padronizando o trace
  P.scaled <- scale.trace(P, target.trace)
  
  d.scaled <- apply(beta, 2, responseDiff,
                    P1 = P.scaled,
                    P2 = Pavg)
  
  results_d <- rbind(results_d,
                     data.frame(
                       Species = sp,
                       Mean_d_raw = mean(d.raw),
                       SD_d_raw = sd(d.raw),
                       Mean_d_scaled = mean(d.scaled),
                       SD_d_scaled = sd(d.scaled)
                     ))
}

saveRDS(results_d, "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Divergence.RDS")
