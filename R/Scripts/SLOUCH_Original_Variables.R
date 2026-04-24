setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(slouch)
library(ape)

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))

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

sem_size_d = vector()
# loop sobre espécies
for(sp in names(by_trait)) {
  macho <- by_trait[[sp]]$Machos       # vetor com 39 traços
  femea <- by_trait[[sp]]$Fêmeas       # vetor com 39 traços
  
  # calcular dimorfismo normalizado
  #dimorfismo_list[[sp]] <- (macho - femea) / ((geomean(macho) + geomean(femea)) / 2)
  isom = rep(0.160128154, 39)
  isom = isom / sqrt(sum(isom^2))
  
  no_size_m = macho - sum(macho * isom) * isom
  no_size_f = femea - sum(femea * isom) * isom
  
  d_shape = no_size_m - no_size_f
  sem_size_d[sp] = norma(d_shape)
}

#dimorfismo_mat <- do.call(rbind, dimorfismo_list)
#rownames(dimorfismo_mat) <- names(dimorfismo_list)
sem_size_d <- sem_size_d[tree$tip.label]

#####################################################
# A PARTIR DAQUI
# Catarrhini
align_c = align[1:38,]
evolvas_c = evolvas[1:38,]
mds_c = mds[1:38,]
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_c$matings.especies))
sem_size_d_c = sem_size_d[1:38]

#
df = data.frame(mds_c$especies, mds_c$dados.PARVORDER, mds_c$dados.SOCIAL_ORGANIZATION,
                mds_c$dados.MATING_SYSTEM, mds_c$dados.AGGRESSION, evolvas_c$Dimorfism, 
                evolvas_c$Integration, evolvas_c$Size, sem_size_d_c)

df <- df[match(tree_c$tip.label, df$mds_c.especies), ]

# Modelos
# BM
df_bm <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(sem_size_d_c, center = TRUE, scale = TRUE)
  response = sem_size_d_c
)

sygma = seq(0.3, 1.1, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_c,
  response = df_bm$response,
  species = df_bm$species,
  sigma2_y_values = sygma
)
plot(fit_bm)

#save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_c,
  response = df_bm$response,
  species =  df_bm$species,
  a_values = alpha_vals
)
plot(fit_ou)

#save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU.RData")

# OU SOCIAL ORG
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(df$align_c.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_c,
  social = df$mds_c.dados.SOCIAL_ORGANIZATION
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_social <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$social,
  a_values = alpha_vals
)
plot(fit_ou_social)

#save(fit_ou_social, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1.RData")

# OU MATING SYSTEM
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(df$align_c.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_c,
  matings = df$mds_c.dados.MATING_SYSTEM
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_matings <- slouch.fit(
  phy = tree_c,
  #response = df_ou$response,
  response = sem_size_d_c,
  species =  df_ou$species,
  random.cov = df_ou$matings
)
plot(fit_ou_matings)

#save(fit_ou_matings, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS2.RData")

# OU AGGRESSION
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(df$align_c.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_c,
  aggression = df$mds_c.dados.AGGRESSION
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_aggression <- slouch.fit(
  phy = tree_c,
  #response = df_ou$response,
  response = sem_size_d_c,
  species =  df_ou$species,
  random.cov = df_ou$aggression
)
plot(fit_ou_aggression)

#save(fit_ou_aggression, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS2.RData")

# OU Integration
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(df$align_c.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_c,
  #integration = scale(df$evolvas_c.Integration, center = TRUE, scale = TRUE)
  integration = df$evolvas_c.Integration
)

alpha_vals <- seq(0.01, 0.4, length.out = 20)
fit_ou_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$integration,
  a_values = alpha_vals
)
plot(fit_ou_integration)
summary(fit_ou_integration)

save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Integration.RData")

# OU SOCIAL + MATINGS
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(df$align_c.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_c,
  social = df$mds_c.dados.SOCIAL_ORGANIZATION,
  matings = df$mds_c.dados.MATING_SYSTEM
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_social_matings <- slouch.fit(
  phy = tree_c,
  #response = df_ou$response,
  response = sem_size_d_c,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "matings")],
  a_values = alpha_vals
)
plot(fit_ou_social_matings)

#save(fit_ou_social_matings, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2.RData")

# OU SOCIAL + MATINGS + AGGRESSION
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(df$align_c.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_c,
  social = df$mds_c.dados.SOCIAL_ORGANIZATION,
  matings = df$mds_c.dados.MATING_SYSTEM,
  aggression = df$mds_c.dados.AGGRESSION
)

alpha_vals <- seq(0.01, 0.2, length.out = 20)
fit_ou_social_matings_aggression <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "matings", "aggression")],
  a_values = alpha_vals
)
plot(fit_ou_social_matings_aggression)

#save(fit_ou_social_matings_aggression, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU SOCIAL + MATINGS + AGGRESSION + INTEGRATION
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  #response = scale(df$align_c.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_c,
  social = df$mds_c.dados.SOCIAL_ORGANIZATION,
  matings = df$mds_c.dados.MATING_SYSTEM,
  aggression = df$mds_c.dados.AGGRESSION,
  integration = df$evolvas_c.Integration
)

alpha_vals <- seq(0.01, 0.3, length.out = 20)
fit_ou_social_matings_aggression_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "matings", "aggression", "integration")],
  a_values = alpha_vals
)
plot(fit_ou_social_matings_aggression_integration)

#save(fit_ou_social_matings_aggression_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

fits = list(
  BM = fit_bm,
  OU = fit_ou,
  OU_SOCIAL = fit_ou_social,
  OU_MATINGS = fit_ou_matings,
  OU_AGGRESSION = fit_ou_aggression,
  OU_INTEGRATION = fit_ou_integration,
  OU_SOCIAL_MATINGS = fit_ou_social_matings,
  OU_SOCIAL_MATINGS_AGGRESSION = fit_ou_social_matings_aggression,
  OU_SOCIAL_MATINGS_AGGRESSION_INTEGRATION = fit_ou_social_matings_aggression_integration
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

# Platyrrhini
align_p = align[39:62,]
evolvas_p = evolvas[39:62,]
mds_p = mds[39:62,]
tree_p = drop.tip(tree, setdiff(tree$tip.label, align_p$matings.especies))
sem_size_d_p = sem_size_d[39:62]

#
df = data.frame(mds_p$especies, mds_p$dados.PARVORDER, mds_p$dados.SOCIAL_ORGANIZATION,
                mds_p$dados.MATING_SYSTEM, mds_p$dados.AGGRESSION, evolvas_p$Dimorfism, 
                evolvas_p$Integration, evolvas_p$Size, sem_size_d_p, mds_p$dados.PROP_MALES_FEMALES)


df <- df[match(tree_p$tip.label, df$mds_p.especies), ]

# Modelos
# BM
df_bm <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE)
  response = sem_size_d_p
)

sygma = seq(0.3, 0.9, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_p,
  response = df_bm$response,
  species = df_bm$species,
  sigma2_y_values = sygma
)
plot(fit_bm)

#save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.04, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_p,
  response = df_bm$response,
  species =  df_bm$species,
  a_values = alpha_vals
)
plot(fit_ou)

#save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU.RData")

# OU SOCIAL (SIZE COMO DUMMY)
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  social = df$mds_p.dados.SOCIAL_ORGANIZATION,
  dummy = df$evolvas_p.Size
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_social <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "dummy")],
  a_values = alpha_vals
)
plot(fit_ou_social)

#save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1.RData")

# OU MATINGS
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  matings = df$mds_p.dados.MATING_SYSTEM
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_matings <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$matings
)
plot(fit_ou_matings)

#save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS2.RData")

# OU PROP M/F
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  prop = scale(df$mds_p.dados.PROP_MALES_FEMALES, center = TRUE, scale = TRUE)
)

alpha_vals <- seq(0.01, 0.07, length.out = 20)
fit_ou_prop <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$prop,
  a_values = alpha_vals
)
plot(fit_ou_prop)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Integration.RData")

# OU AGGRESSION SIZE COMO DUMMY
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  aggression = df$mds_p.dados.AGGRESSION,
  dummy = df$evolvas_p.Size
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_aggression <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("aggression", "dummy")],
  a_values = alpha_vals
)
plot(fit_ou_aggression)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Integration.RData")

# OU INTEGRATION NAO FUNFOU
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  integration = scale(df$evolvas_p.Integration, center = TRUE, scale = TRUE),
  dummy = df$evolvas_p.Size
)

alpha_vals <- seq(0.001, 0.1, length.out = 20)
fit_ou_integration <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("integration", "dummy")],
  a_values = alpha_vals
)
plot(fit_ou_integration)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Integration.RData")

# OU SOCIAL + MATINGS SIZE COMO DUMMY
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  social = df$mds_p.dados.SOCIAL_ORGANIZATION,
  matings = df$mds_p.dados.MATING_SYSTEM,
  dummy = df$evolvas_p.Size
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_social_matings <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "matings", "dummy")],
  a_values = alpha_vals
)
plot(fit_ou_social_matings)

#save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2.RData")

# OU SOCIAL + MATINGS + PROP M/F
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  social = df$mds_p.dados.SOCIAL_ORGANIZATION,
  matings = df$mds_p.dados.MATING_SYSTEM,
  prop = scale(df$mds_p.dados.PROP_MALES_FEMALES, center = TRUE, scale = TRUE)
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_social_matings_prop <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "matings", "prop")],
  a_values = alpha_vals
)
plot(fit_ou_social_matings_prop)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU SOCIAL + MATINGS + PROP M/F + AGGRESSION
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  social = df$mds_p.dados.SOCIAL_ORGANIZATION,
  matings = df$mds_p.dados.MATING_SYSTEM,
  prop = scale(df$mds_p.dados.PROP_MALES_FEMALES, center = TRUE, scale = TRUE),
  aggression = df$mds_p.dados.AGGRESSION
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_social_matings_prop_aggression <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "matings", "prop", "aggression")],
  a_values = alpha_vals
)
plot(fit_ou_social_matings_prop_aggression)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU SOCIAL + MATINGS + PROP M/F + AGGRESSION + INTEGRATION
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  #response = scale(df$align_p.dimor, center = TRUE, scale = TRUE),
  response = sem_size_d_p,
  social = df$mds_p.dados.SOCIAL_ORGANIZATION,
  matings = df$mds_p.dados.MATING_SYSTEM,
  prop = scale(df$mds_p.dados.PROP_MALES_FEMALES, center = TRUE, scale = TRUE),
  aggression = df$mds_p.dados.AGGRESSION,
  integration = scale(df$evolvas_p.Integration, center = TRUE, scale = TRUE)
)

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_social_matings_prop_aggression_integration <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("social", "matings", "prop", "aggression", "integration")],
  a_values = alpha_vals
)
plot(fit_ou_social_matings_prop_aggression_integration)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

fits = list(
  BM = fit_bm,
  OU = fit_ou,
  OU_SOCIAL = fit_ou_social, # SIZE COMO DUMMY
  OU_MATINGS = fit_ou_matings,
  OU_PROP = fit_ou_prop,
  OU_AGGRESSION = fit_ou_aggression, # SIZE COMO DUMMY
  OU_INTEGRATION = fit_ou_integration,
  OU_SOCIAL_MATINGS = fit_ou_social_matings, # SIZE COMO DUMMY
  OU_SOCIAL_MATINGS_PROP = fit_ou_social_matings_prop,
  OU_SOCIAL_MATINGS_PROP_AGGRESSION = fit_ou_social_matings_prop_aggression,
  OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION = fit_ou_social_matings_prop_aggression_integration
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fits.RData")
