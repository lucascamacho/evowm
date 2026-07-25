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

# read all species VCV matrices
setwd("~/Dropbox/Doc/Data/vcv/")
temp <- list.files(pattern = "*.csv")
vcv <- lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv) <- gsub(".csv", replacement = "", temp)
vcv$Lagothrix_lagothricha <- as.data.frame(lapply(vcv$Lagothrix_lagothricha, 
                                                  function(x) 
                                                    as.numeric(trimws(x))))
vcv$Cacajao_calvus <- read.csv("~/Dropbox/Doc/Data/p_vcv_gabriel/Cacajao_calvus.csv", 
                              header = FALSE, sep = ";", 
                              dec = ",")
vcv$Cacajao_calvus <- as.data.frame(lapply(vcv$Cacajao_calvus, 
                                           function(x) 
                                             as.numeric(trimws(x))))
# read phylogeny
mds <- readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
filename <- "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree <- read.nexus(filename)
tree$tip.label[tree$tip.label == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"
species <- mds$especies
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

# get ancestral VCV eigenvectors
all_cov_matrices <- PhyloW(tree, vcv)
ancestral <- getMRCA(tree, tree$tip.label)
pcs <- eigen(all_cov_matrices$`66`)$vectors

# read G matrix
g_m <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Saguinus.rds")$G_cov

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

alpha_vals <- seq(0.01, 0.1, length.out = 20)
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

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_eu_geo <- slouch.fit(
  phy = tree,
  response = df_ou$eudist,
  species =  df_ou$species,
  direct.cov = df_ou$geodist,
  a_values = alpha_vals
)
plot(fit_ou_eu_geo)















