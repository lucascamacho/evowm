setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(evolqg)
library(ape)
library(ggplot2)
library(slouch)
library(dplyr)


# --------------------------------------------------
# Function to remove the first eigenvalue
# --------------------------------------------------

zero_first_eigen <- function(mat) {
  
  eig <- eigen(mat, symmetric = TRUE)
  
  D <- diag(eig$values)
  V <- eig$vectors
  
  D[1, 1] <- max(eig$values) * 1e-8
  
  V %*% D %*% t(V)
}


# --------------------------------------------------
# Geodesic distance
# --------------------------------------------------

geodesic_distance <- function(P1, P2) {
  
  C <- solve(P1, P2)
  
  lambda <- eigen(
    C,
    only.values = TRUE
  )$values
  
  sqrt(sum(log(lambda)^2))
}


# --------------------------------------------------
# Euclidean distance
# --------------------------------------------------

euclidean_distance <- function(P1, P2) {
  
  sqrt(sum((P1 - P2)^2))
}

load(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData"
)

load("~/Desktop/Primatrees.RData")

clades <- read.csv("~/Desktop/taxonomy.csv")

genera <- names(P_matrices$face)

tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, genera)
)

calculate_distances <- function(P_list, tree) {
  
  # Keep only genera present in both the P matrices and tree
  common <- intersect(names(P_list), tree$tip.label)
  
  tree_sub <- drop.tip(
    tree,
    setdiff(tree$tip.label, common)
  )
  
  P_list <- P_list[common]
  
  # Ancestral P matrix
  all_cov_matrices <- PhyloW(tree_sub, P_list)
  
  ancestral <- getMRCA(tree_sub, tree_sub$tip.label)
  
  anc_P <- all_cov_matrices[[as.character(ancestral)]]
  
  # Remove size from ancestral P
  anc_P <- zero_first_eigen(anc_P)
  
  # Remove size from all genus P matrices
  P_list <- lapply(P_list, as.matrix)
  P_list <- lapply(P_list, zero_first_eigen)
  
  # Calculate distances
  results <- data.frame(
    Species = names(P_list),
    GeodesicDistance = vapply(
      P_list,
      geodesic_distance,
      numeric(1),
      P2 = anc_P
    ),
    EuclideanDistance = vapply(
      P_list,
      euclidean_distance,
      numeric(1),
      P2 = anc_P
    )
  )
  
  results
}

results_face <- calculate_distances(
  P_matrices$face,
  tree
)

results_neuro <- calculate_distances(
  P_matrices$neuro,
  tree
)

results_face$Module <- "Face"
results_neuro$Module <- "Neurocranium"

results_g <- rbind(
  results_face,
  results_neuro
)

results_dist <- results_g %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  )

saveRDS(
  results_dist,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Distances_F:N.RDS"
)

df_face <- results_dist %>%
  filter(Module == "Face") %>%
  mutate(
    species = factor(Species, levels = tree$tip.label),
    eudist = EuclideanDistance,
    geodist = GeodesicDistance
  )

df_neuro <- results_dist %>%
  filter(Module == "Neurocranium") %>%
  mutate(
    species = factor(Species, levels = tree$tip.label),
    eudist = EuclideanDistance,
    geodist = GeodesicDistance
  )

tree_face <- drop.tip(
  tree,
  setdiff(tree$tip.label, df_face$species)
)

tree_neuro <- drop.tip(
  tree,
  setdiff(tree$tip.label, df_neuro$species)
)

df_face <- df_face[
  match(tree_face$tip.label, df_face$species),
]

df_neuro <- df_neuro[
  match(tree_neuro$tip.label, df_neuro$species),
]

alpha_vals <- seq(0.1, 0.9, length.out = 20)

fit_face <- slouch.fit(
  phy = tree_face,
  response = df_face$geodist,
  species = df_face$species,
  direct.cov = df_face$eudist,
  a_values = alpha_vals
)

fit_neuro <- slouch.fit(
  phy = tree_neuro,
  response = df_neuro$geodist,
  species = df_neuro$species,
  direct.cov = df_neuro$eudist,
  a_values = alpha_vals
)

summary(fit_face)
face_intercept <- 18.75046
face_slope <- 99.98495
face_r2 <- 0.343

summary(fit_neuro)
neuro_intercept <- 20.64792
neuro_slope <- 54.61225
neuro_r2 <- 0.0573

# Add clade information
results_dist <- results_g %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  )

# Data for connecting Face and Neurocranium within each genus
connections <- results_dist %>%
  select(Species, Module, EuclideanDistance, GeodesicDistance) %>%
  tidyr::pivot_wider(
    names_from = Module,
    values_from = c(EuclideanDistance, GeodesicDistance)
  )

slouch_lines <- data.frame(
  Module = c("Face", "Neurocranium"),
  intercept = c(face_intercept, neuro_intercept),
  slope = c(face_slope, neuro_slope)
)

face_label <- sprintf(
  "Face: slope = %.2f\nR² = %.3f",
  face_slope,
  face_r2
)

neuro_label <- sprintf(
  "Neurocranium: slope = %.2f\nR² = %.3f",
  neuro_slope,
  neuro_r2
)

# Plot
p2 <- ggplot(
  results_dist,
  aes(
    x = EuclideanDistance,
    y = GeodesicDistance,
    color = CLADE,
    shape = Module
  )
) +
  
  # Connect Face and Neurocranium within each genus
  geom_segment(
    data = connections,
    aes(
      x = EuclideanDistance_Face,
      y = GeodesicDistance_Face,
      xend = EuclideanDistance_Neurocranium,
      yend = GeodesicDistance_Neurocranium
    ),
    inherit.aes = FALSE,
    color = "grey60",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  
  # Points
  geom_point(size = 3) +
  
  # SLOUCH lines
  geom_abline(
    data = slouch_lines,
    aes(
      intercept = intercept,
      slope = slope,
      linetype = Module
    ),
    linewidth = 1.2
  ) +
  
  scale_linetype_manual(
    name = "SLOUCH model",
    values = c(
      Face = "solid",
      Neurocranium = "dashed"
    )
  ) +
  
  labs(
    x = "Euclidean Distance",
    y = "Geodesic Distance",
    color = "Clade",
    shape = "Module"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      size = 12
    )
  )

p2

x_max <- max(results_dist$EuclideanDistance)
y_min <- min(results_dist$GeodesicDistance)
y_max <- max(results_dist$GeodesicDistance)

p2 <- p2 +
  annotate(
    "text",
    x = x_max,
    y = y_min + 0.10 * (y_max - y_min),
    label = face_label,
    hjust = 1,
    vjust = 0,
    size = 4.5
  ) +
  annotate(
    "text",
    x = x_max,
    y = y_min + 0.02 * (y_max - y_min),
    label = neuro_label,
    hjust = 1,
    vjust = 0,
    size = 4.5
  )

p2

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure2_F:N.png",
  plot = p2,
  width = 12,
  height = 7,
  dpi = 300
)
