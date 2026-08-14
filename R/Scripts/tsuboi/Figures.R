#Plots
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(slouch)
library(ggplot2)
library(dplyr)
library(ape)
library(gridExtra)

# read all species VCV matrices
#setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
#temp = list.files(pattern = "*.csv")
#vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
#names(vcv)  = gsub(".csv", replacement= "", temp)
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/vcv.RData")

# read phylogeny
#filename <- "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
#tree <- read.nexus(filename)
load("~/Desktop/Primatrees.RData")
species <- names(vcv)
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

# remove vcv which are not in the phylogeny
#vcv <- vcv[!names(vcv) %in% setdiff(names(vcv), tree$tip.label)]
#species <- names(vcv)

# load data
results_g <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Distances.RDS")
results_d <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Divergence.RDS")
results_e <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Averages.RDS")
results_sd <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Sexual_Dimorphism.RDS")
R <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/R_matrix.RDS")
#medias <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")
medias <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData")
sd <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData")
eigen_df <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/PropVar_Rank.RDS")

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
######## FIGURE 3
# 3.1
# Eigenvalues of R
lambda_R <- eigen(
  R,
  symmetric = TRUE
)$values


# Eigenvalues of P for each genus
lambda_P <- t(sapply(names(vcv), function(sp) {
  
  eigen(
    vcv[[sp]],
    symmetric = TRUE
  )$values
  
}))

# Row names = genus names
rownames(lambda_P) <- names(vcv)


# Prepare plotting data
plot_data <- data.frame(
  P = as.vector(t(lambda_P)),
  R = rep(lambda_R, times = nrow(lambda_P)),
  species = rep(rownames(lambda_P), each = 39),
  eigenvector = rep(1:39, times = nrow(lambda_P))
)


# Plot
p3 <- ggplot(
  plot_data,
  aes(
    x = log10(P),
    y = log10(R),
    group = species,
    color = species
  )
) +
  
  geom_point(
    alpha = 0.5
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    alpha = 0.1
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none"
  ) +
  
  labs(
    x = expression(log[10](lambda[P])),
    y = expression(log[10](lambda[R]))
  )

p3

# 3.2
R_diag <- diag(R)
vcv <- lapply(vcv, as.matrix)
P_diag <- lapply(vcv, diag)

RP_slopes <- sapply(names(vcv), function(sp){
  P_diag_sp <- vcv[[sp]] |> diag()
  model <- lm(log(R_diag) ~ log(P_diag_sp))
  coef(model)[2]
  
})

df <- data.frame(RP_slopes = RP_slopes)

p4 = ggplot(df, aes(x = RP_slopes)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20,
                 fill = "steelblue4",
                 color = "white",
                 linewidth = 0.3,
                 alpha = 0.8) +
  geom_density(color = "firebrick",
               linewidth = 1.2) +
  geom_vline(xintercept = mean(RP_slopes),
             linetype = "dashed",
             color = "black",
             linewidth = 0.8) +
  annotate("text",
           x = mean(RP_slopes),
           y = Inf,
           label = sprintf("Mean = %.2f", mean(RP_slopes)),
           vjust = 1.5,
           hjust = -0.1,
           size = 5) +
  labs(
    x = "Slope of log(R) ~ log(P)",
    y = "Density"
  ) +
  theme_classic(base_size = 14)

p4

# 3.3
RP_results <- sapply(names(vcv), function(sp){
  P_diag <- diag(vcv[[sp]])
  model <- lm(log10(R_diag) ~ log10(P_diag))
  c(
    slope = coef(model)[2],
    R2 = summary(model)$r.squared
  )
})

RP_results <- t(RP_results)
RP_results <- as.data.frame(RP_results)

names(results_g$GeodesicDistance)
names(RP_results)[1] <- "RP_slope"
data_integration <- merge(
  RP_results,
  results_g[, c("Species", "GeodesicDistance")],
  by.x = "row.names",
  by.y = "Species"
)

# rename species column
names(data_integration)[1] <- "Species"

data_integration <- data.frame(
  Species = rownames(RP_results),
  GeodesicDistance = results_g$GeodesicDistance[
    match(rownames(RP_results), results_g$Species)
  ],
  RP_slope = RP_results$RP_slope,
  RP_R2 = RP_results$R2
)

df_ou <- data.frame(
  species = factor(data_integration$Species, levels = tree$tip.label),
  rp_r2 = data_integration$RP_R2,
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)
fit_ou_r2_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$rp_r2,
  species =  df_ou$species,
  direct.cov = df_ou$geodist,
  a_values = alpha_vals
)
summary(fit_ou_r2_geodist)

slope <- -0.03437537
intercept <- 1.103648
r2 <- 0.347

p5 <- ggplot(data_integration,
             aes(
               x = GeodesicDistance,
               y = RP_R2
             )) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept,
    slope = slope,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Geodesic distance",
    y = expression(R^2~"of log"[10]*"(R) ~ vs ~ log"[10]*"(P)")
  )

p5

df_ou <- data.frame(
  species = factor(data_integration$Species, levels = tree$tip.label),
  rp_slope = data_integration$RP_slope,
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)
fit_ou_rp_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$rp_slope,
  species =  df_ou$species,
  direct.cov = df_ou$geodist,
  a_values = alpha_vals
)
summary(fit_ou_rp_geodist)

slope <- -0.03792882
intercept <- 1.260445
r2 <- 0.468

p6 <- ggplot(data_integration,
             aes(
               x = GeodesicDistance,
               y = RP_slope
             )) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept,
    slope = slope,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Geodesic distance",
    y = expression("Slope of log"[10]*"(R) ~ vs ~ log"[10]*"(P)")
  )

p6

# fit all figures together p3 ate p6
p_final <- grid.arrange(p3, p4, p5, p6)

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure3.png",
  plot = p_final,
  width = 12,
  height = 7,
  dpi = 300
)

######## FIGURE 4
vcv <- lapply(vcv, as.matrix)

P_mean <- Reduce(
  "+",
  vcv
) / length(vcv)

eig_P_mean <- eigen(
  P_mean,
  symmetric = TRUE
)

# Common eigenvectors
G <- eig_P_mean$vectors

# Eigenvalues of mean P
lambda_P_mean <- eig_P_mean$values

e_values <- sapply(
  names(vcv),
  function(sp) {
    
    P <- vcv[[sp]]
    
    diag(
      t(G) %*%
        P %*%
        G
    )
    
  }
)

c_values <- sapply(
  names(vcv),
  function(sp) {
    
    P <- vcv[[sp]]
    
    P_inv <- solve(P)
    
    1 / diag(
      t(G) %*%
        P_inv %*%
        G
    )
    
  }
)

e_values <- t(e_values)
c_values <- t(c_values)

rownames(e_values) <- names(vcv)
rownames(c_values) <- names(vcv)

colnames(e_values) <- paste0("PC", 1:39)
colnames(c_values) <- paste0("PC", 1:39)


# Check dimensions
dim(e_values)
dim(c_values)


cor_e <- cor(
  e_values,
  method = "pearson"
)

cor_c <- cor(
  c_values,
  method = "pearson"
)

plot_matrix <- matrix(
  NA,
  nrow = 39,
  ncol = 39
)

# Add names BEFORE converting to table
rownames(plot_matrix) <- paste0("PC", 1:39)
colnames(plot_matrix) <- paste0("PC", 1:39)


# Lower triangle = unconditional variance (e)
plot_matrix[lower.tri(plot_matrix)] <-
  cor_e[lower.tri(cor_e)]

# Upper triangle = conditional variance (c)
plot_matrix[upper.tri(plot_matrix)] <-
  cor_c[upper.tri(cor_c)]

# Diagonal = 1
diag(plot_matrix) <- 1

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
  levels = paste0("PC", 39:1)
)

heatmap_data$X <- factor(
  heatmap_data$X,
  levels = paste0("PC", 1:39)
)

p_4 <- ggplot(
  heatmap_data,
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
    )
  ) +
  
  labs(
    x = "Eigenvectors of mean P",
    y = "Eigenvectors of mean P",
    fill = "Correlation"
  )


p_4

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure4.png",
  plot = p_4,
  width = 12,
  height = 12,
  dpi = 300
)

######## FIGURE 5
# 5.1
df_ou <- data.frame(
  species = factor(results_e$Species, levels = tree$tip.label),
  mean_evolvability = log(results_e$Mean_Evolvability),
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)
fit_ou_e_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species =  df_ou$species,
  direct.cov = df_ou$mean_evolvability,
  a_values = alpha_vals
)
summary(fit_ou_e_geodist)

slope <- -0.2164738
intercept <- 23.52066
r2 <- 6.72e-04

data_evolvability <- data.frame(
  log_Evolvability = log(results_e$Mean_Evolvability),
  GeodesicDistance = results_g$GeodesicDistance
)

p7 <- ggplot(data_evolvability,
             aes(
               x = log_Evolvability,
               y = GeodesicDistance
             )) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept,
    slope = slope,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log(Evolvability)),
    y = "Geodesic distance"
  )

p7

# 5.2
df_ou <- data.frame(
  species = factor(results_e$Species, levels = tree$tip.label),
  mean_conditional_evolvability = log(results_e$Mean_Conditional_Evolvability),
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.4, length.out = 20)
fit_ou_c_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species =  df_ou$species,
  direct.cov = df_ou$mean_conditional_evolvability,
  a_values = alpha_vals
)
summary(fit_ou_c_geodist)

intercept <- 21.93461
slope <- -0.2579201
r2 <- 0.0262

data_cond_evolvability <- data.frame(
  log_Conditional_Evolvability = log(results_e$Mean_Conditional_Evolvability),
  GeodesicDistance = results_g$GeodesicDistance
)

p8 <- ggplot(data_cond_evolvability,
             aes(
               x = log_Conditional_Evolvability,
               y = GeodesicDistance
             )) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept,
    slope = slope,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log("Conditional evolvability")),
    y = "Geodesic distance"
  )

p8

# merge p6 and p7 together
p_final2 <- grid.arrange(p7, p8)

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure5.png",
  plot = p_final2,
  width = 12,
  height = 7,
  dpi = 300
)


######## FIGURE EVOLVABILITY SD
# X.1 Evolvability SD
df_ou <- data.frame(
  species = factor(results_sd$Species, levels = tree$tip.label),
  evolvability_sd = log(results_sd$Evolvability_SD),
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.2, 0.9, length.out = 20)
fit_ou_esd_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$evolvability_sd,
  a_values = alpha_vals
)

summary(fit_ou_esd_geodist)

# Fill these after running SLOUCH
intercept <- 22.57885
slope <- -0.3008721
r2 <- 0.175
  
data_evolvability_sd <- data.frame(
    log_Evolvability_SD = log(results_sd$Evolvability_SD),
    GeodesicDistance = results_g$GeodesicDistance
)

p_sd1 <- ggplot(
  data_evolvability_sd,
  aes(
    x = log_Evolvability_SD,
    y = GeodesicDistance
  )
) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept,
    slope = slope,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log("Evolvability SD")),
    y = "Geodesic distance"
  )

p_sd1


# X.2 Conditional Evolvability SD
df_ou <- data.frame(
  species = factor(results_sd$Species, levels = tree$tip.label),
  conditional_evolvability_sd = log(results_sd$Conditional_Evolvability_SD),
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.5, 1, length.out = 20)
fit_ou_cesd_geodist <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species = df_ou$species,
  direct.cov = df_ou$conditional_evolvability_sd,
  a_values = alpha_vals
)

summary(fit_ou_cesd_geodist)

# Fill these after running SLOUCH
intercept <- 24.67428
slope <- -5.599747e-06
r2 <- 5.08e-11
  
data_cond_evolvability_sd <- data.frame(
    log_Conditional_Evolvability_SD = log(results_sd$Conditional_Evolvability_SD),
    GeodesicDistance = results_g$GeodesicDistance
)

p_sd2 <- ggplot(
  data_cond_evolvability_sd,
  aes(
    x = log_Conditional_Evolvability_SD,
    y = GeodesicDistance
  )
) +
  
  geom_point(size = 3) +
  
  geom_abline(
    intercept = intercept,
    slope = slope,
    color = "firebrick",
    linewidth = 1.2
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("R² = %.3f", r2),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = expression(log("Conditional evolvability SD")),
    y = "Geodesic distance"
  )

p_sd2


# Merge plots
p_final_sd <- grid.arrange(p_sd1, p_sd2)

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure_SD.png",
  plot = p_final_sd,
  width = 12,
  height = 7,
  dpi = 300
)

# Figure ranks and proportion of variation
p9 <- ggplot(eigen_df, aes(x = Rank, y = log(Prop_variance), group = Species)) +
        geom_line(alpha = 0.4) +
        geom_point(alpha = 0.4, size = 1) +
        labs(
        x = "Eigenvector rank",
        y = "Log prop of variance explained"
        ) +
     theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )

p10 <- ggplot(eigen_df, aes(x = Rank, y = Prop_variance, group = Species)) +
         geom_line(alpha = 0.4) +
         geom_point(alpha = 0.4, size = 1) +
       labs(
         x = "Eigenvector rank",
         y = "Raw prop of variance explained"
       ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )

p_final_propvar <- grid.arrange(p9, p10)

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Prop_Var_Exp_Rank.png",
  plot = p_final_propvar,
  width = 12,
  height = 10,
  dpi = 300
)
