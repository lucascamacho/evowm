#Plots
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(slouch)
library(ggplot2)
library(dplyr)
library(ape)
library(gridExtra)
library(ggpubr)

load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData")

load("~/Desktop/Primatrees.RData")
species <- names(P_matrices[[1]])
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

clades <- read.csv("~/Desktop/taxonomy.csv")

# load data
results_g <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Distances_F:N.RDS")
results_d <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Divergence_F:N.RDS")
results_e <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Averages_F:N.RDS")
results_sd <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Sexual_Dimorphism_F:N.RDS")
R <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/R_matrices_face_neuro.RData")
medias <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData")
sd <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData")
eigen_df <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/PropVar_Rank_F:N.RDS")

# get species evolvability means
results_e <- results_e %>%
  dplyr::filter(Eigenvector >= 1, Eigenvector <= 38) %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Mean_Evolvability = mean(Evolvability, na.rm = TRUE),
    Mean_Conditional_Evolvability = mean(Conditional_Evolvability, na.rm = TRUE),
    .groups = "drop"
  )

# organize data
sp_comuns <- Reduce(intersect, list(
  results_d$Species,
  results_e$Species,
  results_g$Species,
  results_sd$Species
))

results_d <- results_d[results_d$Species %in% sp_comuns, ]
results_e <- results_e[results_e$Species %in% sp_comuns, ]
results_g <- results_g[results_g$Species %in% sp_comuns, ]
results_sd <- results_sd[results_sd$Species %in% sp_comuns, ]

######################### FIGURES #############################################
lambda_R_face <- eigen(
  R_face,
  symmetric = TRUE
)$values

lambda_R_neuro <- eigen(
  R_neuro,
  symmetric = TRUE
)$values

lambda_P_face <- t(
  sapply(
    names(P_matrices$face),
    function(sp) {
      
      eigen(
        P_matrices$face[[sp]],
        symmetric = TRUE
      )$values
      
    }
  )
)

lambda_P_neuro <- t(
  sapply(
    names(P_matrices$neuro),
    function(sp) {
      
      eigen(
        P_matrices$neuro[[sp]],
        symmetric = TRUE
      )$values
      
    }
  )
)

rownames(lambda_P_face) <- names(P_matrices$face)
rownames(lambda_P_neuro) <- names(P_matrices$neuro)

#
plot_data_face <- data.frame(
  P = as.vector(t(lambda_P_face)),
  R = rep(
    lambda_R_face,
    times = nrow(lambda_P_face)
  ),
  species = rep(
    rownames(lambda_P_face),
    each = ncol(lambda_P_face)
  ),
  eigenvalue = rep(
    1:ncol(lambda_P_face),
    times = nrow(lambda_P_face)
  ),
  Module = "Face"
)

plot_data_neuro <- data.frame(
  P = as.vector(t(lambda_P_neuro)),
  R = rep(
    lambda_R_neuro,
    times = nrow(lambda_P_neuro)
  ),
  species = rep(
    rownames(lambda_P_neuro),
    each = ncol(lambda_P_neuro)
  ),
  eigenvalue = rep(
    1:ncol(lambda_P_neuro),
    times = nrow(lambda_P_neuro)
  ),
  Module = "Neurocranium"
)

plot_data <- rbind(
  plot_data_face,
  plot_data_neuro
)


# Add clade information
plot_data <- plot_data %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("species" = "GENUS")
  )

p3 <- ggplot(
  plot_data,
  aes(
    x = log10(P),
    y = log10(R),
    group = species,
    color = Module
  )
) +
  
  geom_point(
    alpha = 0.5
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  
  facet_wrap(
    ~ Module
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log[10](lambda[P])),
    y = expression(log[10](lambda[R])),
    color = "Module"
  ) +
  
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
  )


p3 # mudar cor para clades

RP_slopes_face <- sapply(
  split(
    plot_data_face,
    plot_data_face$species
  ),
  function(dat) {
    
    model <- lm(
      log10(R) ~ log10(P),
      data = dat
    )
    
    coef(model)[2]
  }
)


RP_slopes_neuro <- sapply(
  split(
    plot_data_neuro,
    plot_data_neuro$species
  ),
  function(dat) {
    
    model <- lm(
      log10(R) ~ log10(P),
      data = dat
    )
    
    coef(model)[2]
  }
)

df_slopes <- rbind(
  
  data.frame(
    species = names(RP_slopes_face),
    RP_slopes = as.numeric(RP_slopes_face),
    Module = "Face"
  ),
  
  data.frame(
    species = names(RP_slopes_neuro),
    RP_slopes = as.numeric(RP_slopes_neuro),
    Module = "Neurocranium"
  )
)


# Add clades
df_slopes <- df_slopes %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("species" = "GENUS")
  )

p4 <- ggplot(
  df_slopes,
  aes(
    x = RP_slopes,
    fill = Module
  )
) +
  
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 20,
    color = "white",
    linewidth = 0.3,
    alpha = 0.8
  ) +
  
  geom_density(
    linewidth = 1.2
  ) +
  
  geom_vline(
    data = df_slopes %>%
      group_by(Module) %>%
      summarise(
        mean_slope = mean(RP_slopes),
        .groups = "drop"
      ),
    aes(
      xintercept = mean_slope
    ),
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ Module
  ) +
  
  labs(
    x = expression(
      "Slope of " ~
        log[10](lambda[R]) ~
        " vs " ~
        log[10](lambda[P])
    ),
    y = "Density"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
  )


p4

calculate_RP_diagonal <- function(R, P_list) {
  
  R_diag <- diag(R)
  
  results <- lapply(
    names(P_list),
    function(sp) {
      
      P_diag <- diag(P_list[[sp]])
      
      model <- lm(
        log10(R_diag) ~ log10(P_diag)
      )
      
      data.frame(
        Species = sp,
        RP_slope = coef(model)[2],
        RP_R2 = summary(model)$r.squared
      )
    }
  )
  
  bind_rows(results)
}


# ---------------------------------------------------------
# Calculate separately
# ---------------------------------------------------------

RP_results_face <- calculate_RP_diagonal(
  R_face,
  P_matrices$face
)

RP_results_neuro <- calculate_RP_diagonal(
  R_neuro,
  P_matrices$neuro
)


# Add module
RP_results_face$Module <- "Face"
RP_results_neuro$Module <- "Neurocranium"


RP_results <- rbind(
  RP_results_face,
  RP_results_neuro
)


# ---------------------------------------------------------
# Add geodesic distance
# ---------------------------------------------------------

data_integration <- RP_results %>%
  left_join(
    results_g %>%
      select(
        Species,
        Module,
        GeodesicDistance
      ),
    by = c(
      "Species",
      "Module"
    )
  )


# Add clade
data_integration <- data_integration %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  )


# ---------------------------------------------------------
# Check
# ---------------------------------------------------------

head(data_integration)

# =========================================================
# SLOUCH: RP_R2 ~ GEODESIC DISTANCE
# =========================================================

df_face <- data_integration %>%
  filter(Module == "Face") %>%
  mutate(
    species = factor(
      Species,
      levels = tree$tip.label
    ),
    rp_r2 = RP_R2,
    geodist = GeodesicDistance
  )


df_neuro <- data_integration %>%
  filter(Module == "Neurocranium") %>%
  mutate(
    species = factor(
      Species,
      levels = tree$tip.label
    ),
    rp_r2 = RP_R2,
    geodist = GeodesicDistance
  )


# Match tree order
df_face <- df_face[
  match(tree$tip.label, df_face$species),
]

df_neuro <- df_neuro[
  match(tree$tip.label, df_neuro$species),
]


# ---------------------------------------------------------
# Fit SLOUCH
# ---------------------------------------------------------

alpha_vals <- seq(
  0.1,
  0.5,
  length.out = 20
)


fit_face_r2 <- slouch.fit(
  phy = tree,
  response = df_face$rp_r2,
  species = df_face$species,
  direct.cov = df_face$geodist,
  a_values = alpha_vals
)


fit_neuro_r2 <- slouch.fit(
  phy = tree,
  response = df_neuro$rp_r2,
  species = df_neuro$species,
  direct.cov = df_neuro$geodist,
  a_values = alpha_vals
)

summary(fit_face_r2)
summary(fit_neuro_r2)

slope_face_r2 <- -0.001055936
intercept_face_r2 <- 0.1096388
r2_face_r2 <- 1.18e-03

slope_neuro_r2 <- -0.01954289
intercept_neuro_r2 <- 0.7269496
r2_neuro_r2 <- 0.144

# =========================================================
# RP SLOPE VS GEODESIC DISTANCE
# =========================================================

df_face <- data_integration %>%
  filter(Module == "Face") %>%
  mutate(
    species = factor(
      Species,
      levels = tree$tip.label
    ),
    rp_slope = RP_slope,
    geodist = GeodesicDistance
  )


df_neuro <- data_integration %>%
  filter(Module == "Neurocranium") %>%
  mutate(
    species = factor(
      Species,
      levels = tree$tip.label
    ),
    rp_slope = RP_slope,
    geodist = GeodesicDistance
  )


df_face <- df_face[
  match(tree$tip.label, df_face$species),
]

df_neuro <- df_neuro[
  match(tree$tip.label, df_neuro$species),
]


# ---------------------------------------------------------
# SLOUCH
# ---------------------------------------------------------

fit_face_slope <- slouch.fit(
  phy = tree,
  response = df_face$rp_slope,
  species = df_face$species,
  direct.cov = df_face$geodist,
  a_values = alpha_vals
)


fit_neuro_slope <- slouch.fit(
  phy = tree,
  response = df_neuro$rp_slope,
  species = df_neuro$species,
  direct.cov = df_neuro$geodist,
  a_values = alpha_vals
)

summary(fit_face_slope)
summary(fit_neuro_slope)

slope_face_slopes <- 0.004449833
intercept_face_slopes <- -0.005302986
r2_face_slopes <- 0.0347

slope_neuro_slopes <- 0.003041813
intercept_neuro_slopes <- 0.1288325
r2_neuro_slopes <- 0.0212

slouch_lines_r2 <- data.frame(
  Module = c("Face", "Neurocranium"),
  intercept = c(
    intercept_face_r2,
    intercept_neuro_r2
  ),
  slope = c(
    slope_face_r2,
    slope_neuro_r2
  )
)

labels_r2 <- data.frame(
  Module = c(
    "Face",
    "Neurocranium"
  ),
  label = c(
    sprintf(
      "SLOUCH R² = %.3f",
      r2_face_r2
    ),
    sprintf(
      "SLOUCH R² = %.3f",
      r2_neuro_r2
    )
  )
)

p5 <- ggplot(
  data_integration,
  aes(
    x = GeodesicDistance,
    y = RP_R2,
    color = CLADE
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_abline(
    data = slouch_lines_r2,
    aes(
      intercept = intercept,
      slope = slope
    ),
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  facet_wrap(
    ~ Module
  ) +
  
  geom_text(
    data = labels_r2,
    aes(
      x = Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  labs(
    x = "Geodesic distance",
    y = expression(
      R^2 ~
        "of log"[10] * "(R) ~ vs ~ log"[10] * "(P)"
    ),
    color = "Clade"
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
  )

p5


# =========================================================
# FIGURE 3.4
# RP SLOPE VS GEODESIC DISTANCE
# =========================================================

slouch_lines_slope <- data.frame(
  Module = c(
    "Face",
    "Neurocranium"
  ),
  intercept = c(
    intercept_face_slopes,
    intercept_neuro_slopes
  ),
  slope = c(
    slope_face_slopes,
    slope_neuro_slopes
  )
)

labels_slope <- data.frame(
  Module = c(
    "Face",
    "Neurocranium"
  ),
  label = c(
    sprintf(
      "SLOUCH R² = %.3f",
      r2_face_slopes
    ),
    sprintf(
      "SLOUCH R² = %.3f",
      r2_neuro_slopes
    )
  )
)


p6 <- ggplot(
  data_integration,
  aes(
    x = GeodesicDistance,
    y = RP_slope,
    color = CLADE
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_abline(
    data = slouch_lines_slope,
    aes(
      intercept = intercept,
      slope = slope
    ),
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  facet_wrap(
    ~ Module
  ) +
  
  geom_text(
    data = labels_slope,
    aes(
      x = Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  labs(
    x = "Geodesic distance",
    y = expression(
      "Slope of log"[10] * "(R) ~ vs ~ log"[10] * "(P)"
    ),
    color = "Clade"
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "bold"
    )
  )

p6

p_final <- ggarrange(
  p3,
  p4,
  p5,
  p6,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  font.label = list(
    size = 16,
    face = "bold"
  )
)

p_final


ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure3_F:N.png",
  plot = p_final,
  width = 18,
  height = 12,
  dpi = 300
)


######## FIGURE 4
# Correlations among eigenvector-specific unconditional
# and conditional variances of P
# =========================================================
# FUNCTION TO CALCULATE HEATMAP DATA
# =========================================================

make_heatmap_data <- function(P_list) {
  
  # Convert matrices to matrices
  P_list <- lapply(P_list, as.matrix)
  
  # -------------------------------------------------------
  # Mean P
  # -------------------------------------------------------
  
  P_mean <- Reduce(
    "+",
    P_list
  ) / length(P_list)
  
  
  # -------------------------------------------------------
  # Eigenvectors of mean P
  # -------------------------------------------------------
  
  eig_P_mean <- eigen(
    P_mean,
    symmetric = TRUE
  )
  
  G <- eig_P_mean$vectors
  
  
  # -------------------------------------------------------
  # Eigenvector-specific unconditional variances
  # -------------------------------------------------------
  
  e_values <- sapply(
    names(P_list),
    function(sp) {
      
      P <- P_list[[sp]]
      
      diag(
        t(G) %*%
          P %*%
          G
      )
      
    }
  )
  
  
  # -------------------------------------------------------
  # Eigenvector-specific conditional variances
  # -------------------------------------------------------
  
  c_values <- sapply(
    names(P_list),
    function(sp) {
      
      P <- P_list[[sp]]
      
      P_inv <- solve(P)
      
      1 / diag(
        t(G) %*%
          P_inv %*%
          G
      )
      
    }
  )
  
  
  # -------------------------------------------------------
  # Transpose
  # -------------------------------------------------------
  
  e_values <- t(e_values)
  c_values <- t(c_values)
  
  
  # -------------------------------------------------------
  # Correlations among eigenvector-specific variances
  # -------------------------------------------------------
  
  cor_e <- cor(
    e_values,
    method = "pearson"
  )
  
  cor_c <- cor(
    c_values,
    method = "pearson"
  )
  
  
  # -------------------------------------------------------
  # Combine into one matrix
  # -------------------------------------------------------
  
  n_traits <- ncol(e_values)
  
  PC_names <- paste0(
    "PC",
    1:n_traits
  )
  
  plot_matrix <- matrix(
    NA,
    nrow = n_traits,
    ncol = n_traits
  )
  
  rownames(plot_matrix) <- PC_names
  colnames(plot_matrix) <- PC_names
  
  
  # Lower triangle = unconditional variance
  plot_matrix[
    lower.tri(plot_matrix)
  ] <- cor_e[
    lower.tri(cor_e)
  ]
  
  
  # Upper triangle = conditional variance
  plot_matrix[
    upper.tri(plot_matrix)
  ] <- cor_c[
    upper.tri(cor_c)
  ]
  
  
  # Diagonal
  diag(plot_matrix) <- 1
  
  
  # -------------------------------------------------------
  # Convert to plotting data
  # -------------------------------------------------------
  
  heatmap_data <- as.data.frame(
    as.table(plot_matrix)
  )
  
  colnames(heatmap_data) <- c(
    "Y",
    "X",
    "Correlation"
  )
  
  
  heatmap_data$Y <- factor(
    heatmap_data$Y,
    levels = rev(PC_names)
  )
  
  heatmap_data$X <- factor(
    heatmap_data$X,
    levels = PC_names
  )
  
  
  return(heatmap_data)
}


# =========================================================
# FACE
# =========================================================

heatmap_face <- make_heatmap_data(
  P_matrices$face
)

heatmap_face$Module <- "Face"


# =========================================================
# NEUROCRANIUM
# =========================================================

heatmap_neuro <- make_heatmap_data(
  P_matrices$neuro
)

heatmap_neuro$Module <- "Neurocranium"

# =========================================================
# HEATMAP FUNCTION
# =========================================================

make_heatmap <- function(data) {
  
  ggplot(
    data,
    aes(
      x = X,
      y = Y,
      fill = Correlation
    )
  ) +
    
    geom_tile() +
    
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "white"
    ) +
    
    coord_equal() +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      axis.title = element_text(
        size = 14
      )
    ) +
    
    labs(
      x = "Eigenvectors of mean P",
      y = "Eigenvectors of mean P",
      fill = "Correlation"
    )
}

p4_face <- make_heatmap(heatmap_face) +
  ggtitle("Face") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  )

p4_neuro <- make_heatmap(heatmap_neuro) +
  ggtitle("Neurocranium") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  )

p4_final <- grid.arrange(
  p4_face,
  p4_neuro,
  ncol = 2
)


ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure4_F:N.png",
  plot = p4_final,
  width = 12,
  height = 6,
  dpi = 300
)

######## FIGURE 5
# EVOLVABILITY AND CONDITIONAL EVOLVABILITY
# SEPARATELY FOR FACE AND NEUROCRANIUM
# =========================================================
# FACE
# =========================================================

data_evolvability_face <- results_e %>%
  filter(Module == "Face") %>%
  left_join(
    results_g %>%
      filter(Module == "Face") %>%
      select(Species, GeodesicDistance),
    by = "Species"
  )

df_ou <- data.frame(
  species = factor(
    data_evolvability_face$Species,
    levels = tree$tip.label
  ),
  mean_evolvability = log(
    data_evolvability_face$Average_Evolvability
  ),
  geodist = data_evolvability_face$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)

fit_ou_e_face <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$mean_evolvability,
  a_values = alpha_vals
)

summary(fit_ou_e_face)

intercept_face <- 11.10375
slope_face <- -2.303914
r2_face <- 0.046

p7_face <- ggplot(
    data_evolvability_face,
    aes(
      x = log(Average_Evolvability),
      y = GeodesicDistance,
      color = CLADE
    )
  ) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept_face,
    slope = slope_face,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2_face),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log(Evolvability)),
    y = "Geodesic distance",
    color = "Clade"
  )

p7_face

# =========================================================
# FACE — CONDITIONAL EVOLVABILITY
# =========================================================

data_cond_evolvability_face <- results_e %>%
  filter(Module == "Face") %>%
  left_join(
    results_g %>%
      filter(Module == "Face") %>%
      select(Species, GeodesicDistance),
    by = "Species"
  )

df_ou <- data.frame(
  species = factor(
    data_cond_evolvability_face$Species,
    levels = tree$tip.label
  ),
  mean_conditional_evolvability = log(
    data_cond_evolvability_face$Average_Conditional_Evolvability
  ),
  geodist = data_cond_evolvability_face$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.4, length.out = 20)

fit_ou_c_face <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$mean_conditional_evolvability,
  a_values = alpha_vals
)

summary(fit_ou_c_face)

intercept_face <- -6.076768
slope_face <- -1.690055
r2_face <- 0.224
  
p8_face <- ggplot(
  data_cond_evolvability_face,
  aes(
    x = log(Average_Conditional_Evolvability),
    y = GeodesicDistance,
    color = CLADE
  )
) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept_face,
    slope = slope_face,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2_face),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log("Conditional evolvability")),
    y = "Geodesic distance",
    color = "Clade"
  )

p8_face  
  
# =========================================================
# NEUROCRANIUM — EVOLVABILITY
# =========================================================

# =========================================================
# FIGURE 5.1 — NEUROCRANIUM
# Average evolvability vs. geodesic distance
# =========================================================

data_evolvability_neuro <- results_e %>%
  filter(Module == "Neurocranium") %>%
  select(
    Species,
    Average_Evolvability
  ) %>%
  distinct() %>%
  left_join(
    results_g %>%
      filter(Module == "Neurocranium") %>%
      select(
        Species,
        GeodesicDistance
      ),
    by = "Species"
  ) %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  ) %>%
  mutate(
    log_Average_Evolvability = log(Average_Evolvability)
  )

nrow(data_evolvability_neuro)

df_ou <- data.frame(
  species = factor(
    data_evolvability_neuro$Species,
    levels = tree$tip.label
  ),
  average_evolvability = data_evolvability_neuro$log_Average_Evolvability,
  geodist = data_evolvability_neuro$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)

fit_ou_e_neuro_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$average_evolvability,
  a_values = alpha_vals
)

summary(fit_ou_e_neuro_geodist)

intercept_neuro <- 9.673728
slope_neuro <- -2.372176
r2_neuro <- 0.0538

p7_neuro <- ggplot(
  data_evolvability_neuro,
  aes(
    x = log(Average_Evolvability),
    y = GeodesicDistance,
    color = CLADE
  )
) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept_neuro,
    slope = slope_neuro,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2_neuro),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log(Evolvability)),
    y = "Geodesic distance",
    color = "Clade"
  )

p7_neuro

# =========================================================
# NEUROCRANIUM — CONDITIONAL EVOLVABILITY
# =========================================================
data_cond_evolvability_neuro <- results_e %>%
  filter(Module == "Neurocranium") %>%
  select(
    Species,
    Average_Conditional_Evolvability
  ) %>%
  distinct() %>%
  left_join(
    results_g %>%
      filter(Module == "Neurocranium") %>%
      select(
        Species,
        GeodesicDistance
      ),
    by = "Species"
  ) %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  ) %>%
  mutate(
    log_Average_Conditional_Evolvability =
      log(Average_Conditional_Evolvability)
  )

nrow(data_cond_evolvability_neuro)

df_ou <- data.frame(
  species = factor(
    data_cond_evolvability_neuro$Species,
    levels = tree$tip.label
  ),
  mean_conditional_evolvability = log(
    data_cond_evolvability_neuro$Average_Conditional_Evolvability
  ),
  geodist = data_cond_evolvability_neuro$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.4, length.out = 20)

fit_ou_c_neuro <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$mean_conditional_evolvability,
  a_values = alpha_vals
)

summary(fit_ou_c_neuro)

intercept_neuro <- -18.71111
slope_neuro <- -2.390962
r2_neuro <- 0.526

p8_neuro <- ggplot(
  data_cond_evolvability_neuro,
  aes(
    x = log(Average_Conditional_Evolvability),
    y = GeodesicDistance,
    color = CLADE
  )
) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept_neuro,
    slope = slope_neuro,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2_neuro),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log("Conditional evolvability")),
    y = "Geodesic distance",
    color = "Clade"
  )

p8_neuro

######## FIGURE SEXUAL DIMORPHISM EVOLVABILITY
######## FACE + NEUROCRANIUM
# =========================================================
# PREPARE FACE DATA
# =========================================================

data_evolvability_sd_face <- results_sd %>%
  filter(Module == "Face") %>%
  select(
    Species,
    Evolvability_SD
  ) %>%
  left_join(
    results_g %>%
      filter(Module == "Face") %>%
      select(
        Species,
        GeodesicDistance
      ),
    by = "Species"
  ) %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  ) %>%
  mutate(
    log_Evolvability_SD = log(Evolvability_SD)
  )

nrow(data_evolvability_sd_face)


data_cond_evolvability_sd_face <- results_sd %>%
  filter(Module == "Face") %>%
  select(
    Species,
    Conditional_Evolvability_SD
  ) %>%
  left_join(
    results_g %>%
      filter(Module == "Face") %>%
      select(
        Species,
        GeodesicDistance
      ),
    by = "Species"
  ) %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  ) %>%
  mutate(
    log_Conditional_Evolvability_SD =
      log(Conditional_Evolvability_SD)
  )

nrow(data_cond_evolvability_sd_face)


# =========================================================
# FACE TREE
# =========================================================

tree_sd_face <- drop.tip(
  tree,
  setdiff(
    tree$tip.label,
    data_evolvability_sd_face$Species
  )
)

length(tree_sd_face$tip.label)


# =========================================================
# FACE — EVOLVABILITY IN SD DIRECTION
# =========================================================

df_ou <- data.frame(
  species = factor(
    data_evolvability_sd_face$Species,
    levels = tree_sd_face$tip.label
  ),
  evolvability_sd =
    data_evolvability_sd_face$log_Evolvability_SD,
  geodist =
    data_evolvability_sd_face$GeodesicDistance
)

df_ou <- df_ou[
  match(
    tree_sd_face$tip.label,
    df_ou$species
  ),
]

stopifnot(
  all(!is.na(df_ou$species)),
  all(!is.na(df_ou$evolvability_sd)),
  all(!is.na(df_ou$geodist)),
  all(
    as.character(df_ou$species) ==
      tree_sd_face$tip.label
  )
)

alpha_vals <- seq(
  0.2,
  0.9,
  length.out = 20
)

fit_ou_esd_face_geodist <- slouch.fit(
  phy = tree_sd_face,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$evolvability_sd,
  a_values = alpha_vals
)

summary(fit_ou_esd_face_geodist)


# =========================================================
# FACE — CONDITIONAL EVOLVABILITY IN SD DIRECTION
# =========================================================

df_ou <- data.frame(
  species = factor(
    data_cond_evolvability_sd_face$Species,
    levels = tree_sd_face$tip.label
  ),
  conditional_evolvability_sd =
    data_cond_evolvability_sd_face$
    log_Conditional_Evolvability_SD,
  geodist =
    data_cond_evolvability_sd_face$GeodesicDistance
)

df_ou <- df_ou[
  match(
    tree_sd_face$tip.label,
    df_ou$species
  ),
]

stopifnot(
  all(!is.na(df_ou$species)),
  all(!is.na(df_ou$conditional_evolvability_sd)),
  all(!is.na(df_ou$geodist)),
  all(
    as.character(df_ou$species) ==
      tree_sd_face$tip.label
  )
)

alpha_vals <- seq(
  0.5,
  1,
  length.out = 20
)

fit_ou_cesd_face_geodist <- slouch.fit(
  phy = tree_sd_face,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$conditional_evolvability_sd,
  a_values = alpha_vals
)

summary(fit_ou_cesd_face_geodist)

# =========================================================
# PREPARE NEUROCRANIUM DATA
# =========================================================

data_evolvability_sd_neuro <- results_sd %>%
  filter(Module == "Neurocranium") %>%
  select(
    Species,
    Evolvability_SD
  ) %>%
  left_join(
    results_g %>%
      filter(Module == "Neurocranium") %>%
      select(
        Species,
        GeodesicDistance
      ),
    by = "Species"
  ) %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  ) %>%
  mutate(
    log_Evolvability_SD = log(Evolvability_SD)
  )

nrow(data_evolvability_sd_neuro)


data_cond_evolvability_sd_neuro <- results_sd %>%
  filter(Module == "Neurocranium") %>%
  select(
    Species,
    Conditional_Evolvability_SD
  ) %>%
  left_join(
    results_g %>%
      filter(Module == "Neurocranium") %>%
      select(
        Species,
        GeodesicDistance
      ),
    by = "Species"
  ) %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Species" = "GENUS")
  ) %>%
  mutate(
    log_Conditional_Evolvability_SD =
      log(Conditional_Evolvability_SD)
  )

nrow(data_cond_evolvability_sd_neuro)


# =========================================================
# NEUROCRANIUM TREE
# =========================================================

tree_sd_neuro <- drop.tip(
  tree,
  setdiff(
    tree$tip.label,
    data_evolvability_sd_neuro$Species
  )
)

length(tree_sd_neuro$tip.label)


# =========================================================
# NEUROCRANIUM — EVOLVABILITY IN SD DIRECTION
# =========================================================

df_ou <- data.frame(
  species = factor(
    data_evolvability_sd_neuro$Species,
    levels = tree_sd_neuro$tip.label
  ),
  evolvability_sd =
    data_evolvability_sd_neuro$log_Evolvability_SD,
  geodist =
    data_evolvability_sd_neuro$GeodesicDistance
)

df_ou <- df_ou[
  match(
    tree_sd_neuro$tip.label,
    df_ou$species
  ),
]

stopifnot(
  all(!is.na(df_ou$species)),
  all(!is.na(df_ou$evolvability_sd)),
  all(!is.na(df_ou$geodist)),
  all(
    as.character(df_ou$species) ==
      tree_sd_neuro$tip.label
  )
)

alpha_vals <- seq(
  0.2,
  0.9,
  length.out = 20
)

fit_ou_esd_neuro_geodist <- slouch.fit(
  phy = tree_sd_neuro,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$evolvability_sd,
  a_values = alpha_vals
)

summary(fit_ou_esd_neuro_geodist)


# =========================================================
# NEUROCRANIUM — CONDITIONAL EVOLVABILITY IN SD DIRECTION
# =========================================================

df_ou <- data.frame(
  species = factor(
    data_cond_evolvability_sd_neuro$Species,
    levels = tree_sd_neuro$tip.label
  ),
  conditional_evolvability_sd =
    data_cond_evolvability_sd_neuro$
    log_Conditional_Evolvability_SD,
  geodist =
    data_cond_evolvability_sd_neuro$GeodesicDistance
)

df_ou <- df_ou[
  match(
    tree_sd_neuro$tip.label,
    df_ou$species
  ),
]

stopifnot(
  all(!is.na(df_ou$species)),
  all(!is.na(df_ou$conditional_evolvability_sd)),
  all(!is.na(df_ou$geodist)),
  all(
    as.character(df_ou$species) ==
      tree_sd_neuro$tip.label
  )
)

alpha_vals <- seq(
  0.5,
  1,
  length.out = 20
)

fit_ou_cesd_neuro_geodist <- slouch.fit(
  phy = tree_sd_neuro,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$conditional_evolvability_sd,
  a_values = alpha_vals
)

summary(fit_ou_cesd_neuro_geodist)

# =========================================================
# SLOUCH RESULTS
# =========================================================

# ---------------------------------------------------------
# FACE — Evolvability SD
# ---------------------------------------------------------

summary(fit_ou_esd_face_geodist)

intercept_face_e <- 17.89865
slope_face_e <- -0.573077
r2_face_e <- 0.276


# ---------------------------------------------------------
# FACE — Conditional Evolvability SD
# ---------------------------------------------------------

summary(fit_ou_cesd_face_geodist)

intercept_face_c <- 23.41841
slope_face_c <- 0.07098532
r2_face_c <- 0.00518


# ---------------------------------------------------------
# NEUROCRANIUM — Evolvability SD
# ---------------------------------------------------------

summary(fit_ou_esd_neuro_geodist)

intercept_neuro_e <- 20.41719
slope_neuro_e <- -0.15841
r2_neuro_e <- 0.0232


# ---------------------------------------------------------
# NEUROCRANIUM — Conditional Evolvability SD
# ---------------------------------------------------------

summary(fit_ou_cesd_neuro_geodist)

intercept_neuro_c <- 21.51671
slope_neuro_c <- -0.01708003
r2_neuro_c <- 5.73e-04

# =========================================================
# FACE — EVOLVABILITY
# =========================================================

p_sd_face_e <- ggplot(
  data_evolvability_sd_face,
  aes(
    x = log_Evolvability_SD,
    y = GeodesicDistance,
    color = CLADE
  )
) +
  geom_point(size = 3) +
  geom_abline(
    intercept = intercept_face_e,
    slope = slope_face_e,
    color = "firebrick",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "R² = %.3f",
      r2_face_e
    ),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression(log("Evolvability SD")),
    y = "Geodesic distance",
    color = "Clade"
  )


# =========================================================
# FACE — CONDITIONAL EVOLVABILITY
# =========================================================

p_sd_face_c <- ggplot(
  data_cond_evolvability_sd_face,
  aes(
    x = log_Conditional_Evolvability_SD,
    y = GeodesicDistance,
    color = CLADE
  )
) +
  geom_point(size = 3) +
  geom_abline(
    intercept = intercept_face_c,
    slope = slope_face_c,
    color = "firebrick",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "R² = %.3f",
      r2_face_c
    ),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression(
      log("Conditional evolvability SD")
    ),
    y = "Geodesic distance",
    color = "Clade"
  )


# =========================================================
# NEUROCRANIUM — EVOLVABILITY
# =========================================================

p_sd_neuro_e <- ggplot(
  data_evolvability_sd_neuro,
  aes(
    x = log_Evolvability_SD,
    y = GeodesicDistance,
    color = CLADE
  )
) +
  geom_point(size = 3) +
  geom_abline(
    intercept = intercept_neuro_e,
    slope = slope_neuro_e,
    color = "firebrick",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "R² = %.3f",
      r2_neuro_e
    ),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression(log("Evolvability SD")),
    y = "Geodesic distance",
    color = "Clade"
  )


# =========================================================
# NEUROCRANIUM — CONDITIONAL EVOLVABILITY
# =========================================================

p_sd_neuro_c <- ggplot(
  data_cond_evolvability_sd_neuro,
  aes(
    x = log_Conditional_Evolvability_SD,
    y = GeodesicDistance,
    color = CLADE
  )
) +
  geom_point(size = 3) +
  geom_abline(
    intercept = intercept_neuro_c,
    slope = slope_neuro_c,
    color = "firebrick",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf(
      "R² = %.3f",
      r2_neuro_c
    ),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression(
      log("Conditional evolvability SD")
    ),
    y = "Geodesic distance",
    color = "Clade"
  )

# =========================================================
# FINAL FIGURE
# =========================================================
p_sd_face_e <- p_sd_face_e +
  labs(title = "Face") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_sd_face_c <- p_sd_face_c +
  labs(title = "Face") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_sd_neuro_e <- p_sd_neuro_e +
  labs(title = "Neurocranium") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_sd_neuro_c <- p_sd_neuro_c +
  labs(title = "Neurocranium") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p7_face <- p7_face +
  labs(title = "Face") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p8_face <- p8_face +
  labs(title = "Face") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p7_neuro <- p7_neuro +
  labs(title = "Neurocranium") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p8_neuro <- p8_neuro +
  labs(title = "Neurocranium") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


p_final_evolvability <- grid.arrange(
  p7_face,
  p8_face,
  p7_neuro,
  p8_neuro,
  ncol = 2
)

p_final_evolvability

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure_5_F:N.png",
  plot = p_final_evolvability,
  width = 12,
  height = 10,
  dpi = 300
)

p_final_sd <- grid.arrange(
  p_sd_face_e,
  p_sd_face_c,
  p_sd_neuro_e,
  p_sd_neuro_c,
  ncol = 2
)

p_final_sd

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure_SD_F:N.png",
  plot = p_final_sd,
  width = 12,
  height = 10,
  dpi = 300
)

# Figure ranks and proportion of variation
clade_genus <- clades %>%
  select(GENUS, FAMILY, CLADE) %>%
  distinct()

eigen_df_clade <- eigen_df %>%
  left_join(
    clade_genus,
    by = c("Genus" = "GENUS")
  )

# --------------------------------------------------
# FACE — Log proportion of variance
# --------------------------------------------------

eigen_face <- eigen_df_clade %>%
  filter(Module == "Face")

p_face_log <- ggplot(
  eigen_face,
  aes(
    x = Rank,
    y = log(Prop_variance),
    group = Genus,
    color = CLADE
  )
) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.4, size = 1) +
  labs(
    x = "Eigenvector rank",
    y = "Log prop. of variance explained",
    color = "Clade"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text = element_text(size = 14),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(size = 12)
  )


# --------------------------------------------------
# FACE — Raw proportion of variance
# --------------------------------------------------

p_face_raw <- ggplot(
  eigen_face,
  aes(
    x = Rank,
    y = Prop_variance,
    group = Genus,
    color = CLADE
  )
) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.4, size = 1) +
  labs(
    x = "Eigenvector rank",
    y = "Raw prop. of variance explained",
    color = "Clade"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text = element_text(size = 14),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(size = 12)
  )


# --------------------------------------------------
# NEUROCRANIUM — Log proportion of variance
# --------------------------------------------------

eigen_neuro <- eigen_df_clade %>%
  filter(Module == "Neurocranium")

p_neuro_log <- ggplot(
  eigen_neuro,
  aes(
    x = Rank,
    y = log(Prop_variance),
    group = Genus,
    color = CLADE
  )
) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.4, size = 1) +
  labs(
    x = "Eigenvector rank",
    y = "Log prop. of variance explained",
    color = "Clade"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text = element_text(size = 14),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(size = 12)
  )


# --------------------------------------------------
# NEUROCRANIUM — Raw proportion of variance
# --------------------------------------------------

p_neuro_raw <- ggplot(
  eigen_neuro,
  aes(
    x = Rank,
    y = Prop_variance,
    group = Genus,
    color = CLADE
  )
) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.4, size = 1) +
  labs(
    x = "Eigenvector rank",
    y = "Raw prop. of variance explained",
    color = "Clade"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA
    ),
    panel.grid = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text = element_text(size = 14),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(size = 12)
  )


# --------------------------------------------------
# Combine the four plots
# --------------------------------------------------

p_face_log <- p_face_log +
  ggtitle("Face — log")

p_face_raw <- p_face_raw +
  ggtitle("Face — raw")

p_neuro_log <- p_neuro_log +
  ggtitle("Neurocranium — log")

p_neuro_raw <- p_neuro_raw +
  ggtitle("Neurocranium — raw")

p_final_propvar <- grid.arrange(
  p_face_log,
  p_face_raw,
  p_neuro_log,
  p_neuro_raw,
  ncol = 2
)

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Prop_Var_Exp_Rank_F:N.png",
  plot = p_final_propvar,
  width = 18,
  height = 10,
  dpi = 300
)
