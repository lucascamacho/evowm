setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ggplot2)
library(ape)

# read all species VCV matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
names(vcv)  = gsub(".csv", replacement= "", temp)

# read phylogeny
filename <- "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree <- read.nexus(filename)
species <- names(vcv)
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

# remove vcv which are not in the phylogeny
vcv <- vcv[!names(vcv) %in% setdiff(names(vcv), tree$tip.label)]
species <- names(vcv)
vcv <- lapply(vcv, as.matrix)

# Empty list to store results
eigen_results <- list()

# Loop through each species
for (i in seq_along(vcv)) {
  
  # Get covariance matrix
  P <- vcv[[i]]
  
  # remove size
  eig <- eigen(P)
  D <- diag(eig$values)
  V <- eig$vectors
  D2 <- D
  D2[1,1] <- max(eig$values) * 1e-8
  P <- V %*% D2 %*% t(V)
  
  # Eigen decomposition
  eig <- eigen(P, symmetric = TRUE)
  
  # Eigenvalues
  eigenvalues <- eig$values
  
  # Proportion of total variance explained
  prop_variance <- eigenvalues / sum(eigenvalues)
  
  # Store results
  eigen_results[[i]] <- data.frame(
    Species = names(vcv)[i],
    Rank = seq_along(eigenvalues),
    Eigenvalue = eigenvalues,
    Prop_variance = prop_variance
  )
}

# Combine all species
eigen_df <- do.call(rbind, eigen_results)

# save results
saveRDS(eigen_df, "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/PropVar_Rank.RDS")
