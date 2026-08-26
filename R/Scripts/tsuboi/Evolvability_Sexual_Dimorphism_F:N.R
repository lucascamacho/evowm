# =========================================================
# EVOLVABILITY IN THE SEXUAL DIMORPHISM DIRECTION
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

load(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData"
)

load(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData"
)

load(
  "~/Desktop/Primaset_RawData.RData"
  )

colnames(genus_means_df)[3:ncol(genus_means_df)] <- distn
colnames(sexual_dimorphism_df)[3:ncol(sexual_dimorphism_df)] <- distn

# Trait names directly from the P matrices
face_traits <- colnames(P_matrices$face[[1]])
neuro_traits <- colnames(P_matrices$neuro[[1]])

# Check that SD data contain all traits
stopifnot(
  all(face_traits %in% colnames(sexual_dimorphism_df)),
  all(neuro_traits %in% colnames(sexual_dimorphism_df))
)

# =========================================================
# REMOVE SIZE
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
# CALCULATE SEXUAL-DIMORPHISM EVOLVABILITY
# =========================================================

calculate_sd_evolvability <- function(
    P_list,
    tree,
    sexual_dimorphism_df,
    trait_names
) {
  
  # -------------------------------------------------------
  # Keep genera present in both P matrices and tree
  # -------------------------------------------------------
  original_species <- names(P_list)
  
  species <- Reduce(
    intersect,
    list(
      original_species,
      tree$tip.label,
      sexual_dimorphism_df$Genus
    )
  )
  
  P_list <- P_list[species]
  
  cat(
    "Number of genera included:",
    length(species),
    "\n"
  )
  
  cat(
    "Excluded genera:",
    paste(
      setdiff(original_species, species),
      collapse = ", "
    ),
    "\n"
  )
  
  # -------------------------------------------------------
  # Results
  # -------------------------------------------------------
  
  results_sd <- data.frame(
    Species = species,
    Evolvability_SD = NA_real_,
    Conditional_Evolvability_SD = NA_real_,
    Average_Evolvability = NA_real_,
    Average_Conditional_Evolvability = NA_real_
  )
  
  # -------------------------------------------------------
  # Calculate for each genus
  # -------------------------------------------------------
  
  for (i in seq_along(species)) {
    
    sp <- species[i]
    
    # -----------------------------------------------------
    # Genus P matrix
    # -----------------------------------------------------
    
    covar <- as.matrix(
      P_list[[sp]]
    )
    
    # Remove size
    covar <- zero_first_eigen(
      covar
    )
    
    
    # -----------------------------------------------------
    # Sexual dimorphism vector
    # -----------------------------------------------------
    
    sd_row <- sexual_dimorphism_df %>%
      filter(Genus == sp)
    
    
    if (nrow(sd_row) != 1) {
      stop(
        paste(
          "Could not uniquely find sexual dimorphism for:",
          sp
        )
      )
    }
    
    
    # Extract only traits belonging to this module
    sd <- as.numeric(
      sd_row[1, trait_names]
    )
    
    
    # Check dimensions
    if (length(sd) != nrow(covar)) {
      stop(
        paste(
          "Dimension mismatch for",
          sp,
          ": SD vector has",
          length(sd),
          "traits but P has",
          nrow(covar),
          "traits."
        )
      )
    }
    
    
    # -----------------------------------------------------
    # Average evolvability
    # -----------------------------------------------------
    
    e_mean <- evolvabilityMeans(
      covar
    )[1]
    
    c_mean <- evolvabilityMeans(
      covar
    )[4]
    
    
    # -----------------------------------------------------
    # Evolvability in SD direction
    # -----------------------------------------------------
    
    evo <- evolvabilityBeta(
      covar,
      sd
    )
    # -----------------------------------------------------
    # Store results
    # -----------------------------------------------------
    
    results_sd[i, ] <- list(
      sp,
      evo$e,
      evo$c,
      e_mean,
      c_mean
    )
  }
  return(results_sd)
}




# =========================================================
# FACE
# =========================================================

results_sd_face <- calculate_sd_evolvability(
  P_list = P_matrices$face,
  tree = tree,
  sexual_dimorphism_df = sexual_dimorphism_df,
  trait_names = face_traits
)

results_sd_face$Module <- "Face"


# =========================================================
# NEUROCRANIUM
# =========================================================

results_sd_neuro <- calculate_sd_evolvability(
  P_list = P_matrices$neuro,
  tree = tree,
  sexual_dimorphism_df = sexual_dimorphism_df,
  trait_names = neuro_traits
)

results_sd_neuro$Module <- "Neurocranium"

results_sd <- bind_rows(
  results_sd_face,
  results_sd_neuro
)

clades <- read.csv(
  "~/Desktop/taxonomy.csv"
)

results_sd <- results_sd %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  )

saveRDS(
  results_sd,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Sexual_Dimorphism_F:N.RDS"
)
