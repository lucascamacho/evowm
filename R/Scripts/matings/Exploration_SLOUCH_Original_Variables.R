setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

library(mvSLOUCH)
library(ape)
library(dplyr)
library(phytools)
library(gt)
library(ggplot2)
library(tidyverse)

# CARREGANDO TODOS OS DADOS
# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")
#nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
  #nicho$genus,
  align$matings.genus,
  evolvas$genus,
  mds$genus
))

#nicho_f   <- nicho[nicho$genus %in% sp_comuns, ]
align_f   <- align[align$matings.genus %in% sp_comuns, ]
evolvas_f <- evolvas[evolvas$genus %in% sp_comuns, ]
mds_f     <- mds[mds$genus %in% sp_comuns, ]

# read all vcv matrices
setwd("~/Dropbox/Doc/Data/genus_vcv/")
temp = list.files(pattern = "*.txt")
vcv = lapply(temp, read.csv, header = FALSE, dec = ".", sep = "\t")
names(vcv)  = gsub(".csv", replacement = "", temp)

wrong <- which(sapply(vcv, ncol) == 1)

for(i in wrong){
  
  file <- temp[i]
  
  # tenta outras combinações
  mat <- read.table(file, header = FALSE, dec = ".")
  
  vcv[[i]] <- as.matrix(mat)
}
vcv <- lapply(vcv, as.matrix)
names(vcv) <- gsub("\\.txt$", "", temp)

# read and plot phylo tree
commom <- intersect(names(vcv), mds_f$genus)
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
genus <- sub("_.*", "", tree$tip.label)
keep <- genus %in% commom
tree2 <- drop.tip(tree, tree$tip.label[!keep])
tree <- tree2
# extrair generos
genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
unique_genera <- unique(genera)
representatives <- c()
for(g in unique_genera){
  
  spp <- tree$tip.label[genera == g]
  
  if(length(spp) == 1){
    
    representatives[g] <- spp
    
  } else {
    
    # pegar comprimento do ramo terminal
    edges <- match(spp, tree$tip.label)
    terminal_lengths <- tree$edge.length[match(edges, tree$edge[,2])]
    
    # escolher o menor (mais recente)
    representatives[g] <- spp[which.max(terminal_lengths)]
  }
}

# remover todas as outras espécies
tree_genus <- drop.tip(tree, setdiff(tree$tip.label, representatives))

# renomear para os gêneros
tree_genus$tip.label <- names(representatives)
plot(tree_genus)

#
#media_genero <- nicho %>%
#  group_by(genus) %>%
#  summarise(
#    climacv_mean = mean(climacv, na.rm = TRUE),
#    CVn_mean = mean(CVn, na.rm = TRUE)
#  )

#nicho_f <- media_genero[match(tree_genus$tip.label, media_genero$genus), ]
evolvas_f <- evolvas_f[match(tree_genus$tip.label, evolvas_f$genus), ]
align_f <- align_f[match(tree_genus$tip.label, align_f$matings.genus), ]
mds_f <- mds_f[match(tree_genus$tip.label, mds_f$genus), ]

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

# RESULTS
# Catarrhini
cata = load("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/catarrhini_fits.RData")

aicc_catar <- c(
  BM = fits$BM$modfit$AICc,
  SIMPLE_OU = fits$OU$modfit$AICc,
  OU_SIZE = fits$OU_SIZE$modfit$AICc,
  OU_NMDS1 = fits$OU_NMDS1$modfit$AICc,
  OU_NMDS2 = fits$OU_NMDS2$modfit$AICc,
  OU_INTEGRATION = fits$OU_INTEGRATION$modfit$AICc,
  #OU_CLIMA = fits$OU_CLIMA$modfit$AICc,
  #OU_CVN = fits$OU_CVN$modfit$AICc,
  OU_NMDS1_NMDS2 = fits$OU_NMDS1_NMDS2$modfit$AICc,
  OU_NMDS1_NMDS2_INTEGRATION = fits$OU_NMDS1_NMDS2_INTEGRATION$modfit$AICc,
  OU_NMDS1_NMDS2_INTEGRATION_SIZE = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE$modfit$AICc
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA$modfit$AICc,
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN$modfit$AICc
)


tab_catar <- make_aicc_table(aicc_catar)

cat("\n=== Catarrhini ===\n")
print(tab_catar, digits = 3)

# Platyrrhini
platy = load("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/platyrrhini_fits.RData")

aicc_platy = c(BM = fits$BM$modfit$AICc,
               SIMPLE_OU = fits$OU$modfit$AICc,
               OU_SIZE = fits$OU_SIZE$modfit$AICc,
               OU_NMDS1 = fits$OU_NMDS1$modfit$AICc,
               OU_NMDS2 = fits$OU_NMDS2$modfit$AICc,
               OU_INTEGRATION = fits$OU_INTEGRATION$modfit$AICc,
               #OU_CLIMA = fits$OU_CLIMA$modfit$AICc,
               #OU_CVN = fits$OU_CVN$modfit$AICc,
               OU_NMDS1_NMDS2 = fits$OU_NMDS1_NMDS2$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION = fits$OU_NMDS1_NMDS2_INTEGRATION$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION_SIZE = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE$modfit$AICc
               #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA$modfit$AICc,
               #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN$modfit$AICc
) 


tab_platy <- make_aicc_table(aicc_platy)

cat("\n=== Platyrrhini ===\n")
print(tab_platy, digits = 3)

# Haplorrhini
haplo = load("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/haplorrhini_fits.RData")

aicc_haplo = c(BM = fits$BM$modfit$AICc,
               SIMPLE_OU = fits$OU$modfit$AICc,
               OU_SIZE = fits$OU_SIZE$modfit$AICc,
               OU_NMDS1 = fits$OU_NMDS1$modfit$AICc,
               OU_NMDS2 = fits$OU_NMDS2$modfit$AICc,
               OU_INTEGRATION = fits$OU_INTEGRATION$modfit$AICc,
               #OU_CLIMA = fits$OU_CLIMA$modfit$AICc,
               #OU_CVN = fits$OU_CVN$modfit$AICc,
               OU_NMDS1_NMDS2 = fits$OU_NMDS1_NMDS2$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION = fits$OU_NMDS1_NMDS2_INTEGRATION$modfit$AICc,
               OU_NMDS1_NMDS2_INTEGRATION_SIZE = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE$modfit$AICc
               #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA$modfit$AICc,
               #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fits$OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN$modfit$AICc
) 


tab_haplo <- make_aicc_table(aicc_haplo)

cat("\n=== Haaplorrhini ===\n")
print(tab_haplo, digits = 3)

# GT TABLES
#################################################################################
library(dplyr)
library(gt)

ou_df <- data.frame(
  
  Clade = c(
    "Catarrhini",
    "Platyrrhini",
    "Platyrrhini",
    "Platyrrhini",
    "Haplorrhini"
  ),
  
  Model = c(
    "Integration",
    "Integration",
    "NMDS1 + NMDS2 + Integration + Size",
    "NMDS1 + NMDS2 + Integration",
    "NMDS1 + NMDS2 + Integration + Size"
  ),
  
  AICc = c(
    33.100,
    5.820,
    6.430,
    7.710,
    22.200
  ),
  
  deltaAICc = c(
    0.000,
    0.000,
    0.610,
    1.890,
    0.000
  ),
  
  rate_of_adaptation = c(
    0.09052632,
    47.89474,
    0.08515562,
    1.74149e-06,
    0.06050325
  ),
  
  phylogenetic_half_life = c(
    7.656858,
    0.0144723,
    8.13977,
    398019.7,
    11.45636
  )
)

gt_ou <- ou_df %>%
  
  gt(groupname_col = "Clade") %>%
  
  cols_label(
    Model = "Model",
    AICc = "AICc",
    deltaAICc = "ΔAICc",
    rate_of_adaptation = "Rate of adaptation",
    phylogenetic_half_life = "Phylogenetic half-life"
  ) %>%
  
  fmt_number(
    columns = c(AICc, deltaAICc, rate_of_adaptation, phylogenetic_half_life),
    decimals = 3
  ) %>%
  
  tab_header(
    title = md("**OU model selection across primate clades (ΔAICc < 2)**")
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  tab_options(
    table.align = "center"
  )

gt_ou

gtsave(
  gt_ou,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ou_parvorder.png",
  vwidth = 1600,
  vheight = 500,
  zoom = 2
)

# PLOT
#################################################################################
cata = load("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/catarrhini_fit_OU_CVn.RData")

align_c = align_f[1:16,]
evolvas_c = evolvas_f[1:16,]
mds_c = mds_f[1:16,]
nicho_c = nicho_f[1:16,]
tree_c = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_c$matings.genus))

# ------------------------------
# 1️⃣ Scatter plot trait vs predictor com slopes
# ------------------------------
trait_obs <- align_c$normas
predictor <-  nicho_c$CVn_mean

# best model
summary(fit_ou_cvn)

# Extrair slopes do summary
optimal_slope <- -4.447655
optimal_se    <- 0.8114961
evol_slope    <- -4.447429
evol_se       <- 0.8114538
intercept     <- 1.235577

df_plot <- data.frame(trait_obs, predictor)

model <- lm(trait_obs ~ predictor, data = df_plot)
r2 <- summary(model)$r.squared

p1 = ggplot(df_plot, aes(x = predictor, y = trait_obs)) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2, color = "#1b9e77") +
  geom_abline(intercept=intercept, slope = evol_slope, color = "blue", linetype = "solid", linewidth = 1.2) +
  geom_abline(intercept = intercept, slope = optimal_slope, color = "red", linetype = "dashed", linewidth = 1.2) +
  # regressão linear dos dados
  geom_smooth(method = "lm", se = TRUE,
              color = "black", linewidth = 1.2) +
  labs(x = "CVn", y = "Sexual Dimorphism") +
  theme_classic() +
  theme(  
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16)
  )
  p1 +
  annotate("text",
           x = Inf, y = -Inf,
           label = paste0("R² = ", round(r2, 3)),
           hjust = 4.1, vjust = -10,
           size = 5,
           color = "#1b9e77")

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/OU_Scatter_CVn_Dimorphism.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução

# ------------------------------
# ------------------------------
# best model
summary(fit_ou_cvn)

alpha <- 641.2443
var_stationary <- 0.07857967
half_life <- 0.001137658
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

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/OU_Stationary_HalfLife_Dimorphism.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


#################################################################################
platy = load("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/platyrrhini_fits.RData")

align_p = align_f[17:28,]
evolvas_p = evolvas_f[17:28,]
mds_p = mds_f[17:28,]
tree_p = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_p$matings.genus))
nicho_p = nicho_f[17:28,]

# ------------------------------
# 1️⃣ Scatter plot trait vs predictor com slopes
# ------------------------------
trait_obs <- align_p$normas
predictor <-  evolvas_p$Integration

# best model
#summary(fit_ou_integration)

# Extrair slopes do summary
optimal_slope <- -7.437947
optimal_se    <- 1.643562
evol_slope    <- -7.404202
evol_se       <- 1.636105
intercept     <- -0.4717971

df_plot <- data.frame(trait_obs, predictor)

p1 = ggplot(df_plot, aes(x = predictor, y = trait_obs)) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2, color = "#7570b3") +
  geom_abline(intercept=intercept, slope = evol_slope, color = "blue", linetype = "solid", linewidth = 1.2) +
  geom_abline(intercept = intercept, slope = optimal_slope, color = "red", linetype = "dashed", linewidth = 1.2) +
  
  labs(x = "Integration Index", y = "Magnitude of Sexual Dimorphism") +
  theme_classic() +
  theme(  
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 16)
  )

p1

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Platy_OU_Scatter_Integration_Dimorphism.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução

# ------------------------------
# ------------------------------
# best model
summary(fit_ou_integration)

alpha <- 10.03105
var_stationary <- 0.03554922
half_life <- 0.06910019
t_max <- max(node.depth.edgelength(tree_p))

time_df <- tibble(
  time = seq(0, t_max, length.out = 500)
) %>%
  mutate(
    var_t = var_stationary * (1 - exp(-2 * alpha * time))
  )

p1 = ggplot(time_df, aes(x = time, y = var_t)) +
  geom_line(linewidth = 2, color = "#7570b3") +
  
  
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

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Platy_OU_Stationary_HalfLife_Dimorphism.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


#####
library(phytools)
library(caper)

normas <- setNames(align_f$normas, align_f$matings.genus)
phylosig(tree_genus, normas, method = "lambda", test = TRUE)

integration <- setNames(evolvas_f$Integration, evolvas_f$genus)
phylosig(tree_genus, integration, method = "lambda", test = TRUE)

lm_fit <- lm(align_f$normas ~ evolvas_f$Integration)
resids <- residuals(lm_fit)
names(resids) <- align_f$matings.genus

lambda_resid <- phylosig(tree_genus, resids, method="lambda", test=TRUE)
lambda_resid

data = data.frame(align_f$matings.genus, align_f$normas, evolvas_f$Integration)
comp <- comparative.data(tree_genus, data, names.col="align_f.matings.genus")

pgls_fit <- pgls(align_f.normas ~ evolvas_f.Integration, data = comp, lambda="ML")
summary(pgls_fit)

plot(comp$data$evolvas_f.Integration,
     comp$data$align_f.normas)

abline(pgls_fit)

pgls_fit2 <- pgls(
  evolvas_f.Integration ~ align_f.normas,
  data = comp,
  lambda = "ML"
)

summary(pgls_fit2)

AIC(pgls_fit, pgls_fit2)


# OU colapsou = RESPOSTA INSTANTANEA DOS GENEROS AO FEITO DA INTEGRACAO NO DIMORFISMO SEXUAL
# 
# Poucas espécies + preditor sem sinal filogenético pode explicar
#
# 
# PGLS com λ capturou a filogenia corretamente → slope continua forte, resíduos ok.
# Conclusão: PGLS é suficiente, OU não acrescenta nada confiável nesse caso.

