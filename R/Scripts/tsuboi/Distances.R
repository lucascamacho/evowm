# Calculate Euclidian and Geodesic Distances between P and Average P

# get the species P matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(evolqg)
library(ape)

zero_first_eigen <- function(mat) {
  eig <- eigen(mat)
  
  D <- diag(eig$values)
  V <- eig$vectors
  
  D[1, 1] <- max(eig$values) * 1e-8
  
  V %*% D %*% t(V)
}

# geodesic distance between average P and P
geodesic_distance <- function(P1, P2){
  
  C <- solve(P1, P2)
  
  lambda <- eigen(C, only.values = TRUE)$values
  
  sqrt(sum(log(lambda)^2))
}

# euclidean distance between average P and P
euclidean_distance <- function(P1, P2){
  
  sqrt(sum((P1 - P2)^2))
}

# read all species VCV matrices
# setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
# temp = list.files(pattern = "*.csv")
# vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
# names(vcv)  = gsub(".csv", replacement= "", temp)
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/vcv.RData")

# read phylogeny
# filename <- "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
# tree <- read.nexus(filename)
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

# zero first eigenvalue of anc_vcv
eig <- eigen(anc_vcv)
D <- diag(eig$values)
V <- eig$vectors
D2 <- D
D2[1,1] <- max(eig$values) * 1e-8
anc_vcv <- V %*% D2 %*% t(V)
#anc_vcv <- RemoveSize(anc_vcv)

# all vcv are matrices and remove the first eigenvalue
vcv <- lapply(vcv, as.matrix)
vcv <- lapply(vcv, zero_first_eigen)
#vcv <- lapply(vcv, RemoveSize)

# apply functions
results_g <- data.frame(
  Species = names(vcv),
  GeodesicDistance = vapply(
    vcv,
    geodesic_distance,
    numeric(1),
    P2 = anc_vcv
  ),
  EuclideanDistance = vapply(
    vcv,
    euclidean_distance,
    numeric(1),
    P2 = anc_vcv
  )
)
rownames(results_g) <- NULL

# save
saveRDS(results_g, "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Distances.RDS")
