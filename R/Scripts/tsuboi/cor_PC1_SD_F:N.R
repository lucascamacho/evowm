# Align SD PCmax by Face and Neurocranium
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ape)
library(evolqg)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Load P matrices
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData")

load("~/Desktop/Primaset_RawData.RData")

# Load phylogeny
load("~/Desktop/Primatrees.RData")

# Sexual dimorphism data
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData")

# Clade information
clades <- read.csv("~/Desktop/taxonomy.csv")

# Measure information
measure_info <- data.frame(
  Measure = distn,
  Module = c(
    "Face",
    "Face",
    "Face",
    "Face",
    "Face",
    "Face",
    "Face",
    "Face",
    "Face",
    "Neurocranium",
    "Neurocranium",
    "Face",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Face",
    "Both",
    "Neurocranium",
    "Neurocranium",
    "Face",
    "Face",
    "Face",
    "Face",
    "Face",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Face",
    "Face",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium",
    "Neurocranium"
  )
)

# ---------------------------------------------------------
# Functions
# ---------------------------------------------------------

prod_interno <- function(x, y) {
  sum(x * y)
}

norma <- function(x) {
  sqrt(prod_interno(x, x))
}

corVector <- function(x, y) {
  prod_interno(x, y) / (norma(x) * norma(y))
}

# Remove size from a P matrix
remove_size <- function(P) {
  
  eig <- eigen(P, symmetric = TRUE)
  
  D <- diag(eig$values)
  
  D[1, 1] <- max(eig$values) * 1e-8
  
  eig$vectors %*% D %*% t(eig$vectors)
}

# ---------------------------------------------------------
# Trait indices
# ---------------------------------------------------------

face_indices <- which(
  measure_info$Module %in% c("Face", "Both")
)

neuro_indices <- which(
  measure_info$Module %in% c("Neurocranium", "Both")
)

face_names <- measure_info$Measure[face_indices]
neuro_names <- measure_info$Measure[neuro_indices]

# ---------------------------------------------------------
# Prepare trees and P matrices
# ---------------------------------------------------------

calculate_ancestral_P <- function(P_list) {
  
  # Species represented in the P matrices
  species <- names(P_list)
  
  # Prune tree
  tree_module <- drop.tip(
    tree,
    setdiff(tree$tip.label, species)
  )
  
  # Make sure matrices are matrices
  P_list <- lapply(P_list, as.matrix)
  
  # Calculate ancestral P
  all_cov_matrices <- PhyloW(
    tree_module,
    P_list
  )
  
  ancestral <- getMRCA(
    tree_module,
    tree_module$tip.label
  )
  
  anc_P <- all_cov_matrices[[as.character(ancestral)]]
  
  # Remove size
  anc_P <- remove_size(anc_P)
  
  return(anc_P)
}

# ---------------------------------------------------------
# Calculate ancestral P for each module
# ---------------------------------------------------------

anc_P_face <- calculate_ancestral_P(
  P_matrices$face
)

anc_P_neuro <- calculate_ancestral_P(
  P_matrices$neuro
)

# ---------------------------------------------------------
# PCmax for each ancestral P
# ---------------------------------------------------------

pcmax_face <- eigen(
  anc_P_face,
  symmetric = TRUE
)$vectors[, 1]

pcmax_neuro <- eigen(
  anc_P_neuro,
  symmetric = TRUE
)$vectors[, 1]

# ---------------------------------------------------------
# Sexual dimorphism vectors
# ---------------------------------------------------------

sd_face <- sexual_dimorphism_df[
  ,
  as.character(face_indices),
  drop = FALSE
]

sd_neuro <- sexual_dimorphism_df[
  ,
  as.character(neuro_indices),
  drop = FALSE
]

# ---------------------------------------------------------
# Calculate correlation + norm
# ---------------------------------------------------------

results_face <- data.frame(
  Genus = sexual_dimorphism_df$Genus,
  Vector_Correlation = apply(
    sd_face,
    1,
    function(x) corVector(x, pcmax_face)
  ),
  Norm = apply(
    sd_face,
    1,
    norma
  ),
  Module = "Face"
)

results_neuro <- data.frame(
  Genus = sexual_dimorphism_df$Genus,
  Vector_Correlation = apply(
    sd_neuro,
    1,
    function(x) corVector(x, pcmax_neuro)
  ),
  Norm = apply(
    sd_neuro,
    1,
    norma
  ),
  Module = "Neurocranium"
)

# Combine
results_sd_pcmax <- rbind(
  results_face,
  results_neuro
)

# ---------------------------------------------------------
# Add clade information
# ---------------------------------------------------------

results_sd_pcmax <- results_sd_pcmax %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Genus" = "GENUS")
  )

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------
# Create connections between Face and Neurocranium for each genus
connections <- results_sd_pcmax %>%
  select(
    Genus,
    Module,
    Vector_Correlation,
    Norm
  ) %>%
  mutate(
    x = atanh(Vector_Correlation)
  ) %>%
  select(Genus, Module, x, Norm) %>%
  tidyr::pivot_wider(
    names_from = Module,
    values_from = c(x, Norm)
  )


p_sd_pcmax <- ggplot(
  results_sd_pcmax,
  aes(
    x = atanh(Vector_Correlation),
    y = Norm,
    color = CLADE,
    shape = Module
  )
) +
  
  # Connect Face and Neurocranium within each genus
  geom_segment(
    data = connections,
    aes(
      x = x_Face,
      y = Norm_Face,
      xend = x_Neurocranium,
      yend = Norm_Neurocranium
    ),
    inherit.aes = FALSE,
    color = "grey60",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 3
  ) +
  
  geom_text_repel(
    aes(label = Genus),
    size = 3,
    show.legend = FALSE
  ) +
  
  labs(
    x = "Z transformed vector correlation of log SD with ancestral PCmax (size removed)",
    y = "Sexual dimorphism vector norm",
    color = "Clade",
    shape = "Module"
  ) +
  
  theme_classic()

p_sd_pcmax

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/CorVector_PCmax_SD_F:N.png",
  plot = p_sd_pcmax,
  width = 12,
  height = 7,
  dpi = 300
)
