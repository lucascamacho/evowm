# =========================================================
# CALCULATE DIVERGENCE BETWEEN P MATRICES
# SEPARATELY FOR FACE AND NEUROCRANIUM
# =========================================================

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ape)
library(evolqg)
library(dplyr)

# =========================================================
# FUNCTIONS
# =========================================================

# ---------------------------------------------------------
# Response difference between two P matrices
# ---------------------------------------------------------

responseDiff <- function(P1, P2, beta) {
  
  D <- P1 - P2
  
  sqrt(
    as.numeric(
      t(beta) %*% D %*% D %*% beta
    )
  )
}


# ---------------------------------------------------------
# Scale trace of P matrix
# ---------------------------------------------------------

scale.trace <- function(P, target.trace) {
  
  P * target.trace / sum(diag(P))
}


# ---------------------------------------------------------
# Remove first eigenvalue
# ---------------------------------------------------------

zero_first_eigen <- function(P) {
  
  eig <- eigen(
    P,
    symmetric = TRUE
  )
  
  D <- diag(eig$values)
  V <- eig$vectors
  
  D[1, 1] <- max(eig$values) * 1e-8
  
  V %*% D %*% t(V)
}


# =========================================================
# LOAD DATA
# =========================================================

load(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData"
)

load("~/Desktop/Primatrees.RData")

clades <- read.csv(
  "~/Desktop/taxonomy.csv"
)


# =========================================================
# FUNCTION TO CALCULATE DIVERGENCE
# =========================================================

calculate_divergence <- function(P_list, tree, nrep = 1000) {
  
  
  # -------------------------------------------------------
  # Keep genera present in both P matrices and tree
  # -------------------------------------------------------
  
  species <- intersect(
    names(P_list),
    tree$tip.label
  )
  
  tree_sub <- drop.tip(
    tree,
    setdiff(tree$tip.label, species)
  )
  
  P_list <- P_list[species]
  
  
  # -------------------------------------------------------
  # Ancestral P matrix
  # -------------------------------------------------------
  
  all_cov_matrices <- PhyloW(
    tree_sub,
    P_list
  )
  
  ancestral <- getMRCA(
    tree_sub,
    tree_sub$tip.label
  )
  
  anc_P <- all_cov_matrices[[as.character(ancestral)]]
  
  
  # -------------------------------------------------------
  # Remove size from ancestral P
  # -------------------------------------------------------
  
  anc_P <- zero_first_eigen(
    anc_P
  )
  
  
  # -------------------------------------------------------
  # Random beta vectors
  # -------------------------------------------------------
  
  ntraits <- nrow(anc_P)
  
  beta <- matrix(
    rnorm(ntraits * nrep),
    ntraits,
    nrep
  )
  
  beta <- apply(
    beta,
    2,
    function(x) {
      x / sqrt(sum(x^2))
    }
  )
  
  
  # -------------------------------------------------------
  # Target trace
  # -------------------------------------------------------
  
  target.trace <- sum(
    diag(anc_P)
  )
  
  
  # -------------------------------------------------------
  # Calculate divergence for each genus
  # -------------------------------------------------------
  
  results_d <- lapply(
    names(P_list),
    function(sp) {
      
      P <- as.matrix(
        P_list[[sp]]
      )
      
      # Remove size from genus P
      P <- zero_first_eigen(
        P
      )
      
      
      # ---------------------------------------------------
      # Raw divergence
      # ---------------------------------------------------
      
      d.raw <- apply(
        beta,
        2,
        responseDiff,
        P1 = P,
        P2 = anc_P
      )
      
      
      # ---------------------------------------------------
      # Trace-scaled divergence
      # ---------------------------------------------------
      
      P.scaled <- scale.trace(
        P,
        target.trace
      )
      
      d.scaled <- apply(
        beta,
        2,
        responseDiff,
        P1 = P.scaled,
        P2 = anc_P
      )
      
      
      # ---------------------------------------------------
      # Results
      # ---------------------------------------------------
      
      data.frame(
        Species = sp,
        Mean_d_raw = mean(d.raw),
        SD_d_raw = sd(d.raw),
        Mean_d_scaled = mean(d.scaled),
        SD_d_scaled = sd(d.scaled)
      )
      
    }
  )
  
  
  bind_rows(results_d)
}


# =========================================================
# FACE
# =========================================================

results_face <- calculate_divergence(
  P_matrices$face,
  tree,
  nrep = 1000
)

results_face$Module <- "Face"


# =========================================================
# NEUROCRANIUM
# =========================================================

results_neuro <- calculate_divergence(
  P_matrices$neuro,
  tree,
  nrep = 1000
)

results_neuro$Module <- "Neurocranium"


# =========================================================
# COMBINE
# =========================================================

results_d <- bind_rows(
  results_face,
  results_neuro
)


# Add clade information
results_d <- results_d %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  )


# =========================================================
# SAVE
# =========================================================

saveRDS(
  results_d,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Divergence_F:N.RDS"
)
