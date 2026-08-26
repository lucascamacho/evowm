# =========================================================
# EVOLVABILITY AND CONDITIONAL EVOLVABILITY
# SEPARATELY FOR FACE AND NEUROCRANIUM
# =========================================================

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(evolvability)
library(dplyr)
library(ape)
library(evolqg)


# =========================================================
# LOAD DATA
# =========================================================

load(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData"
)

load("~/Desktop/Primatrees.RData")


# =========================================================
# FUNCTION TO REMOVE SIZE
# =========================================================

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
# FUNCTION TO CALCULATE EVOLVABILITY
# =========================================================

calculate_evolvability <- function(P_list, tree) {
  
  
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
  # Reconstruct ancestral P
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
  # Eigenvectors of ancestral P
  # -------------------------------------------------------
  
  pcs <- eigen(
    anc_P,
    symmetric = TRUE
  )$vectors
  
  
  # -------------------------------------------------------
  # Number of calculations
  # -------------------------------------------------------
  
  n <- length(species) * ncol(pcs)
  
  
  # -------------------------------------------------------
  # Results object
  # -------------------------------------------------------
  
  results_e <- data.frame(
    Species = character(n),
    Eigenvector = integer(n),
    Evolvability = numeric(n),
    Conditional_Evolvability = numeric(n),
    Average_Evolvability = numeric(n),
    Average_Conditional_Evolvability = numeric(n)
  )
  
  
  # -------------------------------------------------------
  # Calculate e and c
  # -------------------------------------------------------
  
  k <- 1
  
  for (i in seq_along(species)) {
    
    covar <- as.matrix(
      P_list[[species[i]]]
    )
    
    # Remove size from genus P
    covar <- zero_first_eigen(
      covar
    )
    
    
    # Average evolvability
    e_mean <- evolvabilityMeans(
      covar
    )[1]
    
    c_mean <- evolvabilityMeans(
      covar
    )[4]
    
    
    # Evolvability in ancestral P axes
    for (j in seq_len(ncol(pcs))) {
      
      evo <- evolvabilityBeta(
        covar,
        pcs[, j]
      )
      
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
  
  
  return(results_e)
}


# =========================================================
# FACE
# =========================================================

results_e_face <- calculate_evolvability(
  P_matrices$face,
  tree
)

results_e_face$Module <- "Face"


# =========================================================
# NEUROCRANIUM
# =========================================================

results_e_neuro <- calculate_evolvability(
  P_matrices$neuro,
  tree
)

results_e_neuro$Module <- "Neurocranium"


# =========================================================
# COMBINE
# =========================================================

results_e <- bind_rows(
  results_e_face,
  results_e_neuro
)


# =========================================================
# ADD CLADE
# =========================================================

clades <- read.csv(
  "~/Desktop/taxonomy.csv"
)

results_e <- results_e %>%
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
  results_e,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Averages_F:N.RDS"
)
