setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ggplot2)
library(ape)

# Load P matrices
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData")

# Load phylogeny
load("~/Desktop/Primatrees.RData")


# --------------------------------------------------
# Function to calculate eigenvalues and proportion
# of variance for a set of P matrices
# --------------------------------------------------

calculate_eigen_results <- function(P_list, module) {
  
  eigen_results <- list()
  
  for (i in seq_along(P_list)) {
    
    # Get covariance matrix
    P <- as.matrix(P_list[[i]])
    
    # Eigen decomposition
    eig <- eigen(P, symmetric = TRUE)
    
    # --------------------------------------------------
    # Remove size
    # --------------------------------------------------
    
    D <- diag(eig$values)
    
    # Shrink largest eigenvalue
    D[1, 1] <- max(eig$values) * 1e-8
    
    # Reconstruct matrix
    P_no_size <- eig$vectors %*% D %*% t(eig$vectors)
    
    # --------------------------------------------------
    # Eigen decomposition after removing size
    # --------------------------------------------------
    
    eig <- eigen(P_no_size, symmetric = TRUE)
    
    # Eigenvalues
    eigenvalues <- eig$values
    
    # Proportion of total variance explained
    prop_variance <- eigenvalues / sum(eigenvalues)
    
    # Store results
    eigen_results[[i]] <- data.frame(
      Genus = names(P_list)[i],
      Module = module,
      Rank = seq_along(eigenvalues),
      Eigenvalue = eigenvalues,
      Prop_variance = prop_variance
    )
  }
  
  # Combine all genera
  eigen_df <- do.call(rbind, eigen_results)
  
  return(eigen_df)
}


# --------------------------------------------------
# Calculate for face
# --------------------------------------------------

eigen_face <- calculate_eigen_results(
  P_list = P_matrices$face,
  module = "Face"
)


# --------------------------------------------------
# Calculate for neurocranium
# --------------------------------------------------

eigen_neuro <- calculate_eigen_results(
  P_list = P_matrices$neuro,
  module = "Neurocranium"
)


# --------------------------------------------------
# Combine
# --------------------------------------------------

eigen_df <- rbind(
  eigen_face,
  eigen_neuro
)


# --------------------------------------------------
# Save
# --------------------------------------------------

saveRDS(
  eigen_df,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/PropVar_Rank_F:N.RDS"
)
