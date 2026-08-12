# Calculate R matrix from averages data of skulls

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ape)

# R matrix estimation function
evol.vcv_fast <- function(tree, X) {
  
  # Match species
  X <- X[tree$tip.label, ]
  
  # Remove missing values
  keep <- complete.cases(X)
  X <- X[keep, ]
  
  tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(X)))
  
  # Phylogenetic covariance
  C <- vcv(tree)
  
  # Inverse covariance
  Cinv <- solve(C)
  
  # Center traits
  X <- scale(X, center = TRUE, scale = FALSE)
  
  # Evolutionary covariance
  R <- crossprod(X, Cinv %*% X) / (nrow(X)-1)
  
  return(R)
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

# read species avergaes
#medias <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData")

# calculate means for each species
#species_means <- t(sapply(medias$ByTrait_Averages, function(sp) {
#  (sp$Machos + sp$Fêmeas) / 2
#}))
#
#species_means <- as.matrix(species_means)

# calculate R
C <- vcv(tree, corr = FALSE)

X <- as.matrix(genus_means_df[, -1])
rownames(X) <- genus_means_df$Genus

R <- evol.vcv_fast(tree, X)

saveRDS(R, "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/R_matrix.RDS")
