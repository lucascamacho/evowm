############################################################
# Lagged evolutionary correlation (TNT-like) in R
############################################################

library(ape)
library(phytools)
library(dplyr)

#############################
# 1. INPUTS
#############################
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(ape)
library(dplyr)
library(phytools)

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/ancestrals_averages_PCS_autovalues_primates.RDS")
#nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
  #nicho$`esp3$...1`,
  align$matings.especies,
  evolvas$Species,
  mds$especies
))

#nicho_f   <- nicho[nicho$`esp3$...1` %in% sp_comuns, ]
align_f   <- align[align$matings.especies %in% sp_comuns, ]
evolvas_f <- evolvas[evolvas$Species %in% sp_comuns, ]
mds_f     <- mds[mds$especies %in% sp_comuns, ]

# read all VCV matrices
setwd("~/Dropbox/Doc/Data/vcv/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement = "", temp)
vcv$Lagothrix_lagothricha <- as.data.frame(lapply(vcv$Lagothrix_lagothricha, function(x) as.numeric(trimws(x))))
vcv$Cacajao_calvus = read.csv("~/Dropbox/Doc/Data/p_vcv_gabriel/Cacajao_calvus.csv", header = FALSE, sep = ";", dec = ",")
vcv$Cacajao_calvus <- as.data.frame(lapply(vcv$Cacajao_calvus, function(x) as.numeric(trimws(x))))

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
tree$tip.label[tree$tip.label == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"
species = mds_f$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

#
# media_species <- nicho %>%
#   dplyr::group_by(.data$`esp3$...1`) %>%
#   dplyr::summarise(
#     climacv_mean = mean(.data$climacv, na.rm = TRUE),
#     CVn_mean = mean(.data$CVn, na.rm = TRUE)
#   )

#nicho_f <- media_species[match(tree$tip.label, media_species$`esp3$...1`), ]
evolvas_f <- evolvas_f[match(tree$tip.label, evolvas_f$Species), ]
align_f <- align_f[match(tree$tip.label, align_f$matings.especies), ]
mds_f <- mds_f[match(tree$tip.label, mds_f$especies), ]

################################################################################
# HAPLORRHINI
df <- data.frame(
  species = mds_f$especies,
  parvorder = mds_f$dados.PARVORDER,
  social = mds_f$dados.SOCIAL_ORGANIZATION,
  mating = mds_f$dados.MATING_SYSTEM,
  prop = mds_f$dados.PROP_MALES_FEMALES,
  dominance = mds_f$dados.DOMINANCE,
  aggression = mds_f$dados.AGGRESSION,
  dimorfism = evolvas_f$Dimorfism,
  integration = evolvas_f$Integration,
  size = evolvas_f$Size,
  nmds1 = mds_f$scores.fit....1.,
  nmds2 = mds_f$scores.fit....2.,
  normas = align_f$normas
  #clima = nicho_f$climacv_mean,
  #cvn = nicho_f$CVn_mean
)

# X <- evolvas_f$Integration
# Y <- align_f$normas
# names(X) <- evolvas_f$Species
# names(Y) <- align_f$matings.especies

# X <- evolvas_f$Size
# Y <- align_f$normas
# names(X) <- evolvas_f$Species
# names(Y) <- align_f$matings.especies
# 
# X <- df$nmds1
# Y <- align_f$normas
# names(X) <- df$species
# names(Y) <- align_f$matings.especies
# 
X <- df$nmds2
Y <- align_f$normas
names(X) <- df$species
names(Y) <- align_f$matings.especies

################################################################################
# CATARRHINI
align_c = align_f[1:16,]
evolvas_c = evolvas_f[1:16,]
mds_c = mds_f[1:16,]
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_c$matings.especies))
nicho_c = nicho_f[1:16,]

df <- data.frame(
  species = mds_c$especies,
  parvorder = mds_c$dados.PARVORDER,
  social = mds_c$dados.SOCIAL_ORGANIZATION,
  mating = mds_c$dados.MATING_SYSTEM,
  prop = mds_c$dados.PROP_MALES_FEMALES,
  dominance = mds_c$dados.DOMINANCE,
  aggression = mds_c$dados.AGGRESSION,
  dimorfism = evolvas_c$Dimorfism,
  integration = evolvas_c$Integration,
  size = evolvas_c$Size,
  nmds1 = mds_c$scores.fit....1.,
  nmds2 = mds_c$scores.fit....2.,
  normas = align_c$normas,
  clima = nicho_c$climacv_mean,
  cvn = nicho_c$CVn_mean
)

df <- df[match(tree_c$tip.label, df$species), ]

tree <- tree_c
X <- align_c$normas
Y <- evolvas_c$Integration
names(X) <- align_c$matings.especies
names(Y) <- evolvas_c$Species

################################################################################
# PLATYRRHINI
align_p = align_f[17:28,]
evolvas_p = evolvas_f[17:28,]
mds_p = mds_f[17:28,]
tree_p = drop.tip(tree, setdiff(tree$tip.label, align_p$matings.especies))
nicho_p = nicho_f[17:28,]

#
df <- data.frame(
  species = mds_p$especies,
  parvorder = mds_p$dados.PARVORDER,
  social = mds_p$dados.SOCIAL_ORGANIZATION,
  mating = mds_p$dados.MATING_SYSTEM,
  prop = mds_p$dados.PROP_MALES_FEMALES,
  dominance = mds_p$dados.DOMINANCE,
  aggression = mds_p$dados.AGGRESSION,
  dimorfism = evolvas_p$Dimorfism,
  integration = evolvas_p$Integration,
  size = evolvas_p$Size,
  nmds1 = mds_p$scores.fit....1.,
  nmds2 = mds_p$scores.fit....2.,
  normas = align_p$normas,
  clima = nicho_p$climacv_mean,
  cvn = nicho_p$CVn_mean
)

df <- df[match(tree_p$tip.label, df$species), ]

tree <- tree_p
X <- align_p$normas
Y <- evolvas_p$Integration
names(X) <- align_p$matings.especies
names(Y) <- evolvas_p$Species

###############################################################################
## DELCOR
###############################################################################
anc_X <- fastAnc(tree, X)
anc_Y <- fastAnc(tree, Y)

# Combine tip + internal node states
states_X <- c(X, anc_X)
states_Y <- c(Y, anc_Y)

states_X <- states_X[tree$tip.label]
states_Y <- states_Y[tree$tip.label]

Ntip <- length(tree$tip.label)
Nnode <- tree$Nnode

#############################
# 3. COMPUTE BRANCH INCREMENTS
#############################
get_increments <- function(tree, states){
  parent <- as.character(tree$edge[,1])
  child  <- as.character(tree$edge[,2])
  
  states[child] - states[parent]
}

# X
states_vec <- numeric(Ntip + Nnode)
names(states_vec) <- as.character(1:(Ntip + Nnode))
states_vec[as.character(1:Ntip)] <- states_X[tree$tip.label]
states_vec[as.character((Ntip + 1):(Ntip + Nnode))] <- anc_X
incX <- get_increments(tree, states_vec)

# Y
states_vec <- numeric(Ntip + Nnode)
names(states_vec) <- as.character(1:(Ntip + Nnode))
states_vec[as.character(1:Ntip)] <- states_Y[tree$tip.label]
states_vec[as.character((Ntip + 1):(Ntip + Nnode))] <- anc_Y
incY <- get_increments(tree, states_vec)

#############################
# 4. FIND LAGGED PAIRS (TNT-like logic)
#############################

find_pairs <- function(tree, incx, incy,
                       radius = 4,
                       minx = 1e-4,
                       miny = 1e-4,
                       delfac = 1){
  
  pairs <- data.frame(x = numeric(),
                      y = numeric(),
                      w = numeric())
  
  for(i in 1:length(incx)){
    
    xchange <- incx[i]
    
    if(!is.na(xchange) && abs(xchange) >= minx){
      
      node <- tree$edge[i,2]
      dist <- 0
      
      repeat{
        
        parent_edge <- which(tree$edge[,2] == node)
        
        if(length(parent_edge) == 0) break
        
        parent <- tree$edge[parent_edge,1]
        dist <- dist + 1
        
        ychange <- incy[parent_edge]
        
        if(!is.na(ychange) && (abs(ychange) >= miny || dist > radius)){
          
          w <- 1 / (1 + dist * delfac)
          
          pairs <- rbind(pairs,
                         data.frame(x = xchange,
                                    y = ychange,
                                    w = w))
          
          break
        }
        
        node <- parent
      }
      
    }
    
  }
  
  return(pairs)
}

#############################
# 5. EXTRACT PAIRS
#############################

pairs <- find_pairs(tree, incX, incY)

# Remove any NA pairs just in case
pairs <- na.omit(pairs)

#############################
# 6. WEIGHTED REGRESSION
#############################

fit <- lm(y ~ x, data = pairs, weights = w)
summary(fit)

# Observed correlation
obs_r <- cor(pairs$x, pairs$y)
obs_r
