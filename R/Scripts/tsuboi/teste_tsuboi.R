# Tsuboi 2025 code
# This is a basic code to test my data on Tsuboi et al 2025
#
# get the species P matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

# load packages
library(evolvability)
library(evolqg)
library(ape)
library(dplyr)
library(slouch)
library(phytools)
library(parallel)

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

# get ancestral VCV eigenvectors
all_cov_matrices <- PhyloW(tree, vcv)
ancestral <- getMRCA(tree, tree$tip.label)
anc_vcv <- all_cov_matrices[[as.character(ancestral)]]

# zero first eigenvalue
eig <- eigen(anc_vcv)
D <- diag(eig$values)
V <- eig$vectors
D2 <- D
D2[1,1] <- min(eig$values[-1]) * 1e-8
anc_vcv <- V %*% D2 %*% t(V)

# get pcs
pcs <- eigen(anc_vcv)$vectors

# read G matrix
#g_m <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Saguinus.rds")#$G_cov

# calculating e and c in orthogonal axis defined by ancestral P
n <- length(species) * ncol(pcs)

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
  D2[1,1] <- min(eig$values[-1]) * 1e-8
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

# calculating the divergence between species
Pavg <- Reduce("+", vcv) / length(vcv)

# scale trace of P matrix
scale.trace <- function(P, target.trace){
  P * target.trace / sum(diag(P))
}

# generate 1000 vectors
ntraits <- nrow(Pavg)
nrep <- 1000
beta <- matrix(rnorm(ntraits * nrep), ntraits, nrep)
beta <- apply(beta, 2, function(x) x / sqrt(sum(x^2)))

# function to calculate diff
responseDiff <- function(P1, P2, beta){
  
  D <- P1 - P2
  
  sqrt(as.numeric(t(beta) %*% D %*% D %*% beta))
}

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

vcv <- lapply(vcv, as.matrix)

results_g <- data.frame(
  Species = names(vcv),
  GeodesicDistance = vapply(
    vcv,
    geodesic_distance,
    numeric(1),
    P2 = Pavg
  ),
  EuclideanDistance = vapply(
    vcv,
    euclidean_distance,
    numeric(1),
    P2 = Pavg
  )
)

rownames(results_g) <- NULL

sp_comuns <- Reduce(intersect, list(
  results_d$Species,
  results_e$Species,
  results_g$Species
))

results_d <- results_d[results_d$Species %in% sp_comuns, ]
results_e <- results_e[results_e$Species %in% sp_comuns, ]
results_g <- results_g[results_g$Species %in% sp_comuns, ]

results_mean <- results_e %>%
  filter(Eigenvector >= 1, Eigenvector <= 39) %>%
  group_by(Species) %>%
  summarise(
    Mean_Evolvability = mean(Evolvability, na.rm = TRUE),
    Mean_Conditional_Evolvability = mean(Conditional_Evolvability, na.rm = TRUE),
    .groups = "drop"
  )

# Estimate R matrix
medias <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")

species_means <- t(sapply(medias$ByTrait_Averages, function(sp) {
  (sp$Machos + sp$Fêmeas) / 2
}))

species_means <- as.matrix(species_means)

R <- evol.vcv(tree, species_means)

# trait_sets <- list(
#   traits_1_5   = 1:5,
#   traits_6_10  = 6:10,
#   traits_11_15 = 11:15,
#   traits_16_20 = 16:20,
#   traits_21_25 = 21:25,
#   traits_26_30 = 26:30,
#   traits_31_35 = 31:35,
#   traits_36_38 = 36:38
# )
# 
# n_cores <- max(1, detectCores() - 1)
# 
# R_list <- mclapply(
#   trait_sets,
#   function(idx) {
#     evol.vcv(tree, species_means[, idx, drop = FALSE])
#   },
#   mc.cores = n_cores
# )

#saveRDS(R_list, "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/R_list.RDS")
saveRDS(R, "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/R.RDS")

# FIGURES
# FIGURE 5
plot(log(results_mean$Mean_Evolvability), results_g$GeodesicDistance)
plot(log(results_mean$Mean_Conditional_Evolvability), results_g$GeodesicDistance)

df_ou <- data.frame(
  species = factor(results_mean$Species, levels = tree$tip.label),
  mean_evolvability = log(results_mean$Mean_Evolvability),
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)
fit_ou_e_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$mean_evolvability,
  species =  df_ou$species,
  direct.cov = df_ou$geodist,
  a_values = alpha_vals
)
plot(fit_ou_e_geodist)

df_ou <- data.frame(
  species = factor(results_mean$Species, levels = tree$tip.label),
  mean_conditional_evolvability = log(results_mean$Mean_Conditional_Evolvability),
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.4, length.out = 20)
fit_ou_c_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$mean_conditional_evolvability,
  species =  df_ou$species,
  direct.cov = df_ou$geodist,
  a_values = alpha_vals
)
plot(fit_ou_c_geodist)

# FIGURE 2
plot(results_g$EuclideanDistance, results_g$GeodesicDistance)

df_ou <- data.frame(
  species = factor(results_mean$Species, levels = tree$tip.label),
  eudist = results_g$EuclideanDistance,
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)
fit_ou_eu_geo <- slouch.fit(
  phy = tree,
  response = df_ou$eudist,
  species =  df_ou$species,
  direct.cov = df_ou$geodist,
  a_values = alpha_vals
)
plot(fit_ou_eu_geo)

# USING RAW TRAITS (MEANS) AS AXIS INSTEAD OF DO TRANSFORMATIONS



