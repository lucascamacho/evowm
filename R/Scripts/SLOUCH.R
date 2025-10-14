setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(slouch)
library(ape)

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

#Load and plot the phylogeny
# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = read.nexus(filename)
species = mds$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

#
species62 = align$matings.especies

by_trait <- medias$ByTrait_Averages

by_trait <- by_trait[names(by_trait) %in% species62]

# inicializar lista ou matriz pra guardar o dimorfismo
dimorfismo_list <- list()

# loop sobre espécies
for(sp in names(by_trait)) {
  macho <- by_trait[[sp]]$Machos       # vetor com 39 traços
  femea <- by_trait[[sp]]$Fêmeas       # vetor com 39 traços
  
  # calcular dimorfismo normalizado
  dimorfismo_list[[sp]] <- (macho - femea) / geomean(macho)
}

dimorfismo_mat <- do.call(rbind, dimorfismo_list)
rownames(dimorfismo_mat) <- names(dimorfismo_list)
dimorfismo_mat <- dimorfismo_mat[tree$tip.label, ]

#####################################################
# A PARTIR DAQUI
# Catarrhini
align_c = align[1:38,]
evolvas_c = evolvas[1:38,]
mds_c = mds[1:38,]
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_c$matings.especies))

#
df = data.frame(mds_c$especies, mds_c$scores.fit....1., mds_c$scores.fit....2., mds_c$dados.PARVORDER,
                align_c$dimor, align_c$normas, align_c$align_1, align_c$align_2,
                align_c$align_3, align_c$align_4, align_c$align_5, align_c$align_6,
                evolvas_c$Dimorfism, evolvas_c$Integration_A, evolvas_c$Size, 
                evolvas_c$Evolvability, evolvas_c$Conditional_Evolvability, 
                evolvas_c$Average_Evolvability)

df <- df[match(tree_c$tip.label, df$mds_c.especies), ]

# Modelos
# BM
df_bm <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor
)

sygma = seq(0.0001, 0.001, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_c,
  response = df_bm$response,
  species = df_bm$species,
  sigma2_y_values = sygma
)

save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_c,
  response = df_bm$response,
  species =  df_bm$species,
  a_values = alpha_vals
)

save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU.RData")

# OU NMDS1
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds1 = df$mds_c.scores.fit....1.
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds1 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds1,
  a_values = alpha_vals
)

save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds2 = df$mds_c.scores.fit....2.
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds2
)

save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS2.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds1 = df$mds_c.scores.fit....1.,
  nmds2 = df$mds_c.scores.fit....2.
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2")],
  a_values = alpha_vals
)

save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Cond_Evolvability
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds1 = df$mds_c.scores.fit....1.,
  nmds2 = df$mds_c.scores.fit....2.,
  cond_evolv = df$evolvas_c.Conditional_Evolvability
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds1_nmds2_cond_evolv <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "cond_evolv")],
  a_values = alpha_vals
)

save(fit_ou_nmds1_nmds2_cond_evolv, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Cond_Evolvability.RData")

# OU NMDS1 + NMDS2 + Cond_Evolvability + Integration
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds1 = df$mds_c.scores.fit....1.,
  nmds2 = df$mds_c.scores.fit....2.,
  cond_evolv = df$evolvas_c.Conditional_Evolvability,
  integration = df$evolvas_c.Integration_A
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds1_nmds2_cond_evolv_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "cond_evolv", "integration")],
  a_values = alpha_vals
)

save(fit_ou_nmds1_nmds2_cond_evolv_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Cond_Evolvability_Integration.RData")

fits = list(
  BM = fit_bm,
  OU = fit_ou,
  OU_NMDS1 = fit_ou_nmds1,
  OU_NMDS2 = fit_ou_nmds2,
  OU_NMDS1_NMDS2 = fit_ou_nmds1_nmds2,
  OU_NMDS1_NMDS2_COND_EVOLV = fit_ou_nmds1_nmds2_cond_evolv,
  OU_NMDS1_NMDS2_COND_EVOLV_INTEGR = fit_ou_nmds1_nmds2_cond_evolv_integration
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

# Platyrrhini
align_p = align[39:62,]
evolvas_p = evolvas[39:62,]
mds_p = mds[39:62,]
tree_p = drop.tip(tree, setdiff(tree$tip.label, align_p$matings.especies))

#
df = data.frame(mds_p$especies, mds_p$scores.fit....1., mds_p$scores.fit....2., mds_p$dados.PARVORDER,
                align_p$dimor, align_p$normas, align_p$align_1, align_p$align_2,
                align_p$align_3, align_p$align_4, align_p$align_5, align_p$align_6,
                evolvas_p$Dimorfism, evolvas_p$Integration_A, evolvas_p$Size, 
                evolvas_p$Evolvability, evolvas_p$Conditional_Evolvability, 
                evolvas_p$Average_Evolvability)

df <- df[match(tree_p$tip.label, df$mds_p.especies), ]

# Modelos
# BM
df_bm <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.dimor
)

sygma = seq(0.0001, 0.001, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_p,
  response = df_bm$response,
  species = df_bm$species,
  sigma2_y_values = sygma
)

save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.001, 0.2, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_p,
  response = df_bm$response,
  species =  df_bm$species,
  a_values = alpha_vals
)

save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU.RData")

# OU NMDS1
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.dimor,
  nmds1 = df$mds_p.scores.fit....1.
)

alpha_vals <- seq(0.00001, 0.00005, length.out = 20)
fit_ou_nmds1 <- slouch.fit(
  phy = tree_p,
  response = scale(df_ou$response),
  species =  df_ou$species,
  random.cov = scale(df_ou$nmds1),
  a_values = 0.00004
  
)
# PAREI AQUI

summary(fit_ou_nmds1)
save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds2 = df$mds_c.scores.fit....2.
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds2
)

save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS2.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds1 = df$mds_c.scores.fit....1.,
  nmds2 = df$mds_c.scores.fit....2.
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2")],
  a_values = alpha_vals
)

save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Cond_Evolvability
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds1 = df$mds_c.scores.fit....1.,
  nmds2 = df$mds_c.scores.fit....2.,
  cond_evolv = df$evolvas_c.Conditional_Evolvability
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds1_nmds2_cond_evolv <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "cond_evolv")],
  a_values = alpha_vals
)

save(fit_ou_nmds1_nmds2_cond_evolv, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Cond_Evolvability.RData")

# OU NMDS1 + NMDS2 + Cond_Evolvability + Integration
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.dimor,
  nmds1 = df$mds_c.scores.fit....1.,
  nmds2 = df$mds_c.scores.fit....2.,
  cond_evolv = df$evolvas_c.Conditional_Evolvability,
  integration = df$evolvas_c.Integration_A
)

alpha_vals <- seq(0.01, 0.5, length.out = 20)
fit_ou_nmds1_nmds2_cond_evolv_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "cond_evolv", "integration")],
  a_values = alpha_vals
)

save(fit_ou_nmds1_nmds2_cond_evolv_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Cond_Evolvability_Integration.RData")

fits = list(
  BM = fit_bm,
  OU = fit_ou,
  OU_NMDS1 = fit_ou_nmds1,
  OU_NMDS2 = fit_ou_nmds2,
  OU_NMDS1_NMDS2 = fit_ou_nmds1_nmds2,
  OU_NMDS1_NMDS2_COND_EVOLV = fit_ou_nmds1_nmds2_cond_evolv,
  OU_NMDS1_NMDS2_COND_EVOLV_INTEGR = fit_ou_nmds1_nmds2_cond_evolv_integration
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")
