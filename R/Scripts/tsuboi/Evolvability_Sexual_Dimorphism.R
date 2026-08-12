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

# read species avergaes
#medias <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData")
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData")

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

n = length(species)

# list for results
results_sd <- data.frame(
  Species = character(n),
  Evolvability_SD = numeric(n),
  Conditional_Evolvability_SD = numeric(n),
  Average_Evolvability = numeric(n),
  Average_Conditional_Evolvability = numeric(n)
)


for(i in seq_along(species)){
  
  covar <- as.matrix(vcv[[species[i]]])
  #medias_sp <- medias$ByTrait_Averages[[species[i]]]
  sd <- as.numeric(
    sexual_dimorphism_df[
      sexual_dimorphism_df$Genus == species[i],
      -(1:3)
    ]
  )
  
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
  
 # sd <- medias_sp$Machos - medias_sp$Fêmeas
  
  evo <- evolvabilityBeta(covar, sd)
  
  results_sd[i, ] <- list(
    species[i],
    evo$e,
    evo$c,
    e_mean,
    c_mean)

}

#
saveRDS(results_sd, "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Sexual_Dimorphism.RDS")

