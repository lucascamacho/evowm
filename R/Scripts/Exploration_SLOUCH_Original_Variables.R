setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(ape)
library(dplyr)
library(phytools)
library(gt)
library(ggplot2)
library(tidyverse)
library(reshape2)

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
tree$tip.label[tree$tip.label == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"
species = mds_f$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

align_filtrado <- align[align$matings.especies %in% mds$especies, ]
evolvas_filtrado <- evolvas[evolvas$Species %in% mds$especies, ]

# RESULTS
# Catarrhini
cata = load("~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

aicc_catar = c(BM = fits$BM$modfit$AICc,
               OU_SIMPLE = fits$OU$modfit$AICc,
               OU_SIZE = fits$OU_SIZE$modfit$AICc,
               OU_NMDS1 = fits$OU_NMDS1$modfit$AICc,
               OU_NMDS2 = fits$OU_NMDS2$modfit$AICc,
               OU_INTEGRATION = fits$OU_INTEGRATION$modfit$AICc,
               #OU_CLIMA = fits$OU_CLIMA$modfit$AICc,
               #OU_CVN = fits$OU_CVN$modfit$AICc,
               OU_NMDS1_NMDS2 = fits$OU_NMDS1_NMDS2$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION = fits$OU_NMDS1_NMDS2_INTEGRATION$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION_SIZE = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE$modfit$AICc
               #OU_OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA$modfit$AICc,
               #OU_OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN$modfit$AICc
)

# -----------------------------
# 2) Função para montar tabela
# -----------------------------
make_aicc_table <- function(aicc_values) {
  data.frame(
    model = names(aicc_values),
    AICc = as.numeric(aicc_values)
  ) %>%
    arrange(AICc) %>%
    mutate(
      deltaAICc = AICc - min(AICc, na.rm = TRUE),
      weight = exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc))
    )
}

tab_catar <- make_aicc_table(aicc_catar)

cat("\n=== Catarrhini ===\n")
print(tab_catar, digits = 3)

#
platy = load("~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fits.RData")

aicc_platy = c(BM = fits$BM$modfit$AICc,
               OU_SIMPLE = fits$OU$modfit$AICc,
               OU_SIZE = fits$OU_SIZE$modfit$AICc,
               OU_NMDS1 = fits$OU_NMDS1$modfit$AICc,
               OU_NMDS2 = fits$OU_NMDS2$modfit$AICc,
               OU_INTEGRATION = fits$OU_INTEGRATION$modfit$AICc,
               #OU_CLIMA = fits$OU_CLIMA$modfit$AICc,
               #OU_CVN = fits$OU_CVN$modfit$AICc,
               OU_NMDS1_NMDS2 = fits$OU_NMDS1_NMDS2$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION = fits$OU_NMDS1_NMDS2_INTEGRATION$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION_SIZE = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE$modfit$AICc
               #OU_OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA$modfit$AICc,
               #OU_OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN$modfit$AICc
)

make_aicc_table <- function(aicc_values) {
  data.frame(
    model = names(aicc_values),
    AICc = as.numeric(aicc_values)
  ) %>%
    arrange(AICc) %>%
    mutate(
      deltaAICc = AICc - min(AICc, na.rm = TRUE),
      weight = exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc))
    )
}

tab_platy <- make_aicc_table(aicc_platy)

cat("\n=== Platyrrhini ===\n")
print(tab_platy, digits = 3)
# MELHOR MODELO E BM, CVN, SIMPLE E NMDS1 (VARIOS MODELOS EXPLICAM)

# Haplo
haplo = load("~/Dropbox/Doc/Code/evowm/R/Outputs/haplorrhini_fits.RData")

aicc_haplo = c(BM = fits$BM$modfit$AICc,
               OU_SIMPLE = fits$OU$modfit$AICc,
               OU_SIZE = fits$OU_SIZE$modfit$AICc,
               OU_NMDS1 = fits$OU_NMDS1$modfit$AICc,
               OU_NMDS2 = fits$OU_NMDS2$modfit$AICc,
               OU_INTEGRATION = fits$OU_INTEGRATION$modfit$AICc,
               OU_CLIMA = fits$OU_CLIMA$modfit$AICc,
               OU_CVN = fits$OU_CVN$modfit$AICc,
               OU_NMDS1_NMDS2 = fits$OU_NMDS1_NMDS2$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION = fits$OU_NMDS1_NMDS2_INTEGRATION$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION_SIZE = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE$modfit$AICc,
               OU_OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA$modfit$AICc,
               OU_OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN$modfit$AICc
)

make_aicc_table <- function(aicc_values) {
  data.frame(
    model = names(aicc_values),
    AICc = as.numeric(aicc_values)
  ) %>%
    arrange(AICc) %>%
    mutate(
      deltaAICc = AICc - min(AICc, na.rm = TRUE),
      weight = exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc))
    )
}

tab_haplo <- make_aicc_table(aicc_haplo)

cat("\n=== Haplorrhini ===\n")
print(tab_haplo, digits = 3)
# MELHOR MODELO E O MAIS COMPLEXO OU_OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN

# PLOTS
#################################################################################
cata = load("~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Size.RData")

# ------------------------------
# 1️⃣ Scatter plot trait vs predictor com slopes
# ------------------------------
trait_obs <- align_filtrado$normas[1:38]                # substitua
predictor <- evolvas_filtrado$Size[1:38]           # seu predictor

summary(fit_ou_size)

# Extrair slopes do summary
optimal_slope <- -42.72985
optimal_se    <- 5.931631
evol_slope    <- -0.05992167
evol_se       <- 0.008318147
intercept     <- 1.213436

df_plot <- data.frame(trait_obs, predictor)

p1 = ggplot(df_plot, aes(x=predictor, y=trait_obs)) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2, color = "#1b9e77") +
  geom_abline(intercept=intercept, slope=evol_slope, color="blue", linetype="solid", linewidth=1.2) +
  geom_abline(intercept=intercept, slope=optimal_slope, color="red", linetype="dashed", linewidth=1.2) +
  labs(x="Species Average Size", y="Sexual Dimorphism") +
  theme_classic() +
  theme(  
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16)
  )

p1

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/OU_Scatter_Integration_Dimorphism.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução

# ------------------------------
# 5
# ------------------------------
summary(fit_ou_size)

alpha <- 9.474186e-05
var_stationary <- 11.70788
half_life <- 7316.166
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_filtrado$matings.especies[1:38]))
t_max <- max(node.depth.edgelength(tree_c))

time_df <- tibble(
  time = seq(0, t_max, length.out = 500)
) %>%
  mutate(
    var_t = var_stationary * (1 - exp(-2 * alpha * time))
  )

p1 = ggplot(time_df, aes(x = time, y = var_t)) +
  geom_line(linewidth = 2, color = "#1b9e77") +
  
  
  # Variância estacionária
  geom_hline(
    yintercept = var_stationary,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  
  # Texto da variância estacionária
  annotate(
    "text",
    x = t_max * 0.6,
    y = var_stationary,
    label = "Stationary variance",
    vjust = -0.6,
    size = 6
  ) +
  
  
  # Half-life
  geom_vline(
    xintercept = half_life,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  
  annotate(
    "text",
    x = half_life,
    y = var_stationary * 0.35,
    label = paste0("Phylogenetic Half-life = ", round(half_life, 2), " Myr"),
    angle = 90,
    vjust = -0.7,
    size = 6
  ) +
  
  # Ponto da variância atual
  #geom_point(
  #  aes(x = t_max, y = 76.77285),
  #  size = 3,
  #  color = "#1b9e77"
  #) +
  
  
  #annotate(
  #  "text",
  #  x = t_max,
  #  y = 76.77285,
  #  label = "Current variance",
  #  hjust = 0.5,
  #  vjust = 3.2,
  #  size = 6
  #) +
  
  
  labs(
    x = "Time since root (Myr)",
    y = "Accumulated trait variance"
  ) +
  coord_cartesian(xlim = c(0, t_max * 1.05)) +
  theme_classic() +
  
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

p1

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/OU_Stationary_HalfLife_Dimorphism.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


#################################################################################
platy = load("~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")
summary(fit_bm)

# Dimorfismo sexual
dimorfismo <- align_filtrado$normas[39:65]  # seu vetor de valores
dimorfismo <- setNames(dimorfismo, align_filtrado$matings.especies[39:65])

tree_p = drop.tip(tree, setdiff(tree$tip.label, align_filtrado$matings.especies[39:65]))
dimorfismo <- dimorfismo[match(tree_p$tip.label, names(dimorfismo))]

n_sims <- 10000
sigma <- sqrt(0.008822727)
intercept <- -0.2768349

p1 = ggplot(data.frame(dimorfismo), aes(x=dimorfismo)) +
  geom_histogram(binwidth=0.5, fill="#7570b3", color="#7570b3", alpha=0.7) +
  geom_vline(xintercept=intercept, color="red", linetype="dashed", linewidth=1.2) +
  labs(x="Sexual Dimorphism", y="Frequency") +
       #title="Distribuição do dimorfismo sexual",
       #subtitle="Linha vermelha = intercepto do modelo BM") +
  theme_classic() +
  theme(  
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

p1

# Simulação BM
sim_traits <- replicate(n_sims, fastBM(tree_p, sig2 = sigma^2, a = intercept))

# Organizando para ggplot
sim_df <- as.data.frame(sim_traits)
sim_df$species <- rownames(sim_df)
sim_melt <- melt(sim_df, id.vars="species", variable.name="simulation", value.name="trait")

# Observados
obs_df <- data.frame(species=tree_p$tip.label, trait=dimorfismo)

# Plot
p1 = ggplot() +
  geom_line(data=sim_melt, aes(x=species, y=trait, group=simulation),
            color="red", alpha=0.2) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2, data=obs_df, aes(x=species, y=trait),
             color="#7570b3") +
  labs(x="Species", y="Sexual Dimorphism") +
       #title="Dimorfismo sob modelo BM",
       #subtitle="Red lines = BM simulations | Blue points = observados") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
        )
p1

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/BM_Scatter_Platyrrhini.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução
# ------------------------------
# 1️⃣ Simulação de traços BM
# ------------------------------
sim_traits <- replicate(n_sims, fastBM(tree_p, sig2=sigma^2, a=intercept))

# Organizar em data.frame para ggplot
sim_df <- as.data.frame(sim_traits)
sim_df$species <- rownames(sim_df)
sim_melt <- melt(sim_df, id.vars="species", variable.name="simulation", value.name="trait")

# Dados observados
obs_df <- data.frame(species = tree_p$tip.label, trait = dimorfismo)

# ------------------------------
# 2️⃣ Boxplot / Violin plot da distribuição simulada
# ------------------------------
p1 = ggplot(sim_melt, aes(x=species, y=trait)) +
  geom_violin(fill="skyblue", alpha=0.2) +
  geom_jitter(data=obs_df, aes(x=species, y=trait), color="#7570b3", size=4, width=0.2, stroke = 2) +
  labs(x="Species", y="Sexual Dimorphism") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 10),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
  )
p1

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/BM_Violin_Platyrrhini.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução
