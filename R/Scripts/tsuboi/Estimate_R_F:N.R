# Calculate R matrices separately for Face and Neurocranium

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ape)

# ---------------------------------------------------------
# R matrix estimation function
# ---------------------------------------------------------

evol.vcv_fast <- function(tree, X) {
  
  # Match species
  X <- X[tree$tip.label, , drop = FALSE]
  
  # Remove species with missing values
  keep <- complete.cases(X)
  X <- X[keep, , drop = FALSE]
  
  # Remove those species from tree
  tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(X)))
  
  # Phylogenetic covariance
  C <- vcv(tree)
  
  # Inverse covariance
  Cinv <- solve(C)
  
  # Center traits
  X <- scale(X, center = TRUE, scale = FALSE)
  
  # Evolutionary covariance
  R <- crossprod(X, Cinv %*% X) / (nrow(X) - 1)
  
  return(R)
}


# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData")

load("~/Desktop/Primatrees.RData")

load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData")

clades <- read.csv("~/Desktop/taxonomy.csv")

load("~/Desktop/Primaset_RawData.RData")

# ---------------------------------------------------------
# Species means
# ---------------------------------------------------------

X <- as.matrix(genus_means_df[, -c(1:2)])
rownames(X) <- genus_means_df$Genus


# ---------------------------------------------------------
# Trait classification
# ---------------------------------------------------------

measure_info <- data.frame(
  Measure = distn,
  Module = c(
    "Face", "Face", "Face", "Face", "Face", "Face",
    "Face", "Face", "Face",
    "Neurocranium", "Neurocranium",
    "Face",
    "Neurocranium", "Neurocranium", "Neurocranium",
    "Neurocranium", "Neurocranium", "Neurocranium",
    "Face",
    "Both",
    "Neurocranium", "Neurocranium",
    "Face", "Face", "Face", "Face", "Face",
    "Neurocranium", "Neurocranium", "Neurocranium",
    "Neurocranium",
    "Face", "Face",
    "Neurocranium", "Neurocranium", "Neurocranium",
    "Neurocranium", "Neurocranium", "Neurocranium"
  )
)


# ---------------------------------------------------------
# Make sure trait names match
# ---------------------------------------------------------

colnames(X) <- measure_info$Measure

# Reorder X according to measure_info
X <- X[, measure_info$Measure, drop = FALSE]


# ---------------------------------------------------------
# Face traits
# ---------------------------------------------------------

face_traits <- measure_info$Measure[
  measure_info$Module %in% c("Face", "Both")
]

X_face <- X[, face_traits, drop = FALSE]


# ---------------------------------------------------------
# Neurocranium traits
# ---------------------------------------------------------

neuro_traits <- measure_info$Measure[
  measure_info$Module %in% c("Neurocranium", "Both")
]

X_neuro <- X[, neuro_traits, drop = FALSE]


# ---------------------------------------------------------
# Phylogeny
# ---------------------------------------------------------

species <- rownames(X)

tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, species)
)


# ---------------------------------------------------------
# Calculate R matrices
# ---------------------------------------------------------

R_face <- evol.vcv_fast(
  tree = tree,
  X = X_face
)

R_neuro <- evol.vcv_fast(
  tree = tree,
  X = X_neuro
)


# ---------------------------------------------------------
# Check dimensions
# ---------------------------------------------------------

dim(R_face)
dim(R_neuro)

colnames(R_face)
colnames(R_neuro)


# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

save(
  R_face,
  R_neuro,
  face_traits,
  neuro_traits,
  file = "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/R_matrices_face_neuro.RData"
)
