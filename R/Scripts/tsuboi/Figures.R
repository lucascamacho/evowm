#Plots
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(slouch)
library(ggplot2)
library(dplyr)
library(ape)

# read all species VCV matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
names(vcv)  = gsub(".csv", replacement= "", temp)

# read phylogeny
filename <- "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree <- read.nexus(filename)
species <- names(vcv)
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

# remove vcv which are not in the phylogeny
vcv <- vcv[!names(vcv) %in% setdiff(names(vcv), tree$tip.label)]
species <- names(vcv)

# load data
results_g <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Distances.RDS")
results_d <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Divergence.RDS")
results_e <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Averages.RDS")
results_sd <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Evolvability_Sexual_Dimorphism.RDS")
R <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/R_matrix.RDS")
medias <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")

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
  results_g$Species
))

results_d <- results_d[results_d$Species %in% sp_comuns, ]
results_e <- results_e[results_e$Species %in% sp_comuns, ]
results_g <- results_g[results_g$Species %in% sp_comuns, ]

######################### FIGURES #############################################
###### FIGURE 2
df_ou <- data.frame(
  species = factor(results_g$Species, levels = tree$tip.label),
  eudist = results_g$EuclideanDistance,
  geodist = results_g$GeodesicDistance
)

df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]

alpha_vals <- seq(0.1, 0.5, length.out = 20)
fit_ou_eu_geo <- slouch.fit(
  phy = tree,
  response = df_ou$geodist,
  species =  df_ou$species,
  direct.cov = df_ou$eudist,
  a_values = alpha_vals
)
summary(fit_ou_eu_geo)

intercept <- 19.78099
slope <- 100.6858
slope_se <- 26.81016
r2 <- 0.162

label <- sprintf("Slope = %.2f \u00B1 %.2f\nR² = %.3f",
                 slope, slope_se, r2)

p2 <- ggplot(results_g, aes(x = EuclideanDistance,
                            y = GeodesicDistance)) +
  geom_point(size = 2.5, alpha = 0.7, color = "steelblue4") +
  geom_abline(intercept = intercept,
              slope = slope,
              color = "firebrick",
              linewidth = 1.2) +
  annotate("text",
           x = Inf, y = -Inf,
           label = label,
           hjust = 1.05, vjust = -0.5,
           size = 5) +
  labs(
    x = "Euclidean Distance",
    y = "Geodesic Distance"
  ) +
  theme_classic(base_size = 14)

p2

######## FIGURE 3
# 3.1
eig_R <- eigen(R, symmetric = TRUE)
lambda_R <- eig_R$values

lambda_P <- t(sapply(medias$Autovalues, function(x){
  unlist(x)
}))

lambda_R <- eigen(R, symmetric = TRUE)$values[1:6]

plot_data <- data.frame(
  P = as.vector(t(lambda_P)),
  R = rep(lambda_R, times = nrow(lambda_P)),
  species = rep(rownames(lambda_P), each = 6),
  eigenvector = rep(1:6, times = nrow(lambda_P))
)


p3 <- ggplot(plot_data, aes(
  x = log10(P),
  y = log10(R)
)) +
  
  geom_point(
    alpha = 0.5
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  
  theme_classic() +
  
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
  
  theme_classic() +
  
  labs(
    x = "Geodesic distance",
    y = "R-P scaling slope"
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
  
  geom_point(size = 3,
             alpha = 0.8,
             color = "steelblue4") +
  
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
    y = expression(R^2~"(R-P relationship)")
  )

p6

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

slope <- -2.413214
intercept <- 9.890194
r2 <- 0.0541

data_evolvability <- data.frame(
  log_Evolvability = log(results_e$Mean_Evolvability),
  GeodesicDistance = results_g$GeodesicDistance
)

p7 <- ggplot(data_evolvability,
             aes(
               x = log_Evolvability,
               y = GeodesicDistance
             )) +
  
  geom_point(size = 3,
             alpha = 0.8,
             color = "steelblue4") +
  
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

slope <- -1.348708
intercept <- 7.514148
r2 <- 0.604

data_cond_evolvability <- data.frame(
  log_Conditional_Evolvability = log(results_e$Mean_Conditional_Evolvability),
  GeodesicDistance = results_g$GeodesicDistance
)

p8 <- ggplot(data_cond_evolvability,
             aes(
               x = log_Conditional_Evolvability,
               y = GeodesicDistance
             )) +
  
  geom_point(size = 3,
             alpha = 0.8,
             color = "steelblue4") +
  
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

######## FIGURE EVOLVABILITY SD


