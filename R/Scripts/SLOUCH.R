# SLOUCH POR ESPECIE
# CATARRHINI, PLATYRRHINI, HAPLORRHINI
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
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/ancestrals_averages_PCS_autovalues_primates.RDS")
#nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
  #nicho$`esp3$...1`,
  align$matings.especies,
  evolvas$Species,
  mds$especies
))

#nicho_f   <- nicho[nicho$`esp3$...1` %in% sp_comuns, ]
align_f   <- align[align$matings.especies %in% sp_comuns, ]
evolvas_f <- evolvas[evolvas$Species %in% sp_comuns, ]
mds_f     <- mds[mds$especies %in% sp_comuns, ]

# read all VCV matrices
setwd("~/Dropbox/Doc/Data/vcv/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement = "", temp)
vcv$Lagothrix_lagothricha <- as.data.frame(lapply(vcv$Lagothrix_lagothricha, function(x) as.numeric(trimws(x))))
vcv$Cacajao_calvus = read.csv("~/Dropbox/Doc/Data/p_vcv_gabriel/Cacajao_calvus.csv", header = FALSE, sep = ";", dec = ",")
vcv$Cacajao_calvus <- as.data.frame(lapply(vcv$Cacajao_calvus, function(x) as.numeric(trimws(x))))

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
tree$tip.label[tree$tip.label == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"
species = mds_f$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

#
#nicho_f <- nicho_f[match(tree$tip.label, nicho_f$`esp3$...1`), ]
evolvas_f <- evolvas_f[match(tree$tip.label, evolvas_f$Species), ]
align_f <- align_f[match(tree$tip.label, align_f$matings.especies), ]
mds_f <- mds_f[match(tree$tip.label, mds_f$especies), ]

################################################################################
# A PARTIR DAQUI
# Catarrhini
align_c = align_f[1:38,]
evolvas_c = evolvas_f[1:38,]
mds_c = mds_f[1:38,]
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_c$matings.especies))
#nicho_c = nicho_f[1:34,]

#
df = data.frame(mds_c$especies, mds_c$dados.PARVORDER, mds_c$dados.SOCIAL_ORGANIZATION,
                mds_c$dados.MATING_SYSTEM, mds_c$dados.AGGRESSION, evolvas_c$Dimorfism, 
                evolvas_c$Integration, evolvas_c$Size, mds_c$scores.fit....1., 
                mds_c$scores.fit....2., align_c$normas)
                #, nicho_c$climacv, nicho_c$CVn)

df <- df[match(tree_c$tip.label, df$mds_c.especies), ]

# Modelos
# BM
df_bm <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = align_c$normas
)
df_bm <- df_bm[match(tree_c$tip.label, df_bm$species), ]

sygma = seq(0.01, 0.03, length.out = 20)
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

# OU Size
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = df$align_c.normas,
  size = scale(df$evolvas_c.Size)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]

fit_ou_size <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$size
)
summary(fit_ou_size)

#save(fit_ou_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_size.RData")

# OU NMDS1
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.scores.fit....1.)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds1 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = NULL,
  a_values = alpha_vals
)
plot(fit_ou_nmds1)

#save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds2 = scale(df$mds_c.scores.fit....2.)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds2,
  a_values = alpha_vals
)
plot(fit_ou_nmds2)

#save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS2.RData")

# OU Integration
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = align_c$normas,
  integration = scale(df$evolvas_c.Integration)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$integration,
  a_values = alpha_vals
)
plot(fit_ou_integration)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Integration.RData")

# # OU Clima
# df_ou <- data.frame(
#   species = factor(df$mds_c.especies, levels = tree_c$tip.label),
#   response = align_c$normas,
#   clima = df$nicho_c.climacv
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]
# 
# alpha_vals <- seq(0.01, 0.1, length.out = 20)
# fit_ou_clima <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou$clima,
#   a_values = alpha_vals
# )
# plot(fit_ou_clima)
# 
# #save(fit_ou_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Clima.RData")
# 
# # OU CVn
# df_ou <- data.frame(
#   species = factor(df$mds_c.especies, levels = tree_c$tip.label),
#   response = align_c$normas,
#   cvn = df$nicho_c.CVn
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]
# 
# fit_ou_cvn <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou$cvn
# )
# summary(fit_ou_cvn)

#save(fit_ou_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_CVn.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.scores.fit....1.),
  nmds2 = scale(df$mds_c.scores.fit....2.)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]

alpha_vals <- seq(0.02, 0.1, length.out = 20)
fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2")],
  a_values = alpha_vals
)
plot(fit_ou_nmds1_nmds2)

#save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Integration
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.scores.fit....1.),
  nmds2 = scale(df$mds_c.scores.fit....2.),
  integration = scale(df$evolvas_c.Integration)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]

alpha_vals <- seq(0.01, 0.3, length.out = 20)
fit_ou_nmds1_nmds2_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "integration")],
  a_values = alpha_vals
)
plot(fit_ou_nmds1_nmds2_integration)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU NMDS1 + NMDS2 + Integration + Size
df_ou <- data.frame(
  species = factor(df$mds_c.especies, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.scores.fit....1.),
  nmds2 = scale(df$mds_c.scores.fit....2.),
  integration = scale(df$evolvas_c.Integration),
  size = scale(df$evolvas_c.Size)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]

fit_ou_nmds1_nmds2_integration_size <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size")]
)
summary(fit_ou_nmds1_nmds2_integration_size)

# # OU NMDS1 + NMDS2 + Integration + Size + Clima
# df_ou <- data.frame(
#   species = factor(df$mds_c.especies, levels = tree_c$tip.label),
#   response = align_c$normas,
#   nmds1 = df$mds_c.scores.fit....1.,
#   nmds2 = df$mds_c.scores.fit....2.,
#   integration = df$evolvas_c.Integration,
#   size = df$evolvas_c.Size,
#   clima = df$nicho_c.climacv
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]
# 
# fit_ou_nmds1_nmds2_integration_size_clima <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima")]
# )
# summary(fit_ou_nmds1_nmds2_integration_size_clima)
# 
# #save(fit_ou_nmds1_nmds2_integration_size_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima.RData")
# 
# # OU NMDS1 + NMDS2 + Integration + Size + Clima + CVn
# df_ou <- data.frame(
#   species = factor(df$mds_c.especies, levels = tree_c$tip.label),
#   response = align_c$normas,
#   nmds1 = df$mds_c.scores.fit....1.,
#   nmds2 = df$mds_c.scores.fit....2.,
#   integration = df$evolvas_c.Integration,
#   size = df$evolvas_c.Size,
#   clima = df$nicho_c.climacv,
#   cvn = df$nicho_c.CVn
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$species), ]
# 
# fit_ou_nmds1_nmds2_integration_size_clima_cvn <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima", "cvn")]
# )
# summary(fit_ou_nmds1_nmds2_integration_size_clima_cvn)
# 
# #save(fit_ou_nmds1_nmds2_integration_size_clima_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima_CVn.RData")

fits = list(
  BM = fit_bm,
  OU = fit_ou,
  OU_SIZE = fit_ou_size,
  OU_NMDS1 = fit_ou_nmds1,
  OU_NMDS2 = fit_ou_nmds2,
  OU_INTEGRATION = fit_ou_integration,
  #OU_CLIMA = fit_ou_clima,
  #OU_CVN = fit_ou_cvn,
  OU_NMDS1_NMDS2 = fit_ou_nmds1_nmds2,
  OU_NMDS1_NMDS2_INTEGRATION = fit_ou_nmds1_nmds2_integration,
  OU_NMDS1_NMDS2_INTEGRATION_SIZE = fit_ou_nmds1_nmds2_integration_size
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fit_ou_nmds1_nmds2_integration_size_clima,
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fit_ou_nmds1_nmds2_integration_size_clima_cvn
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

#################################################################################
# Platyrrhini
align_p = align_f[39:65,]
evolvas_p = evolvas_f[39:65,]
mds_p = mds_f[39:65,]
tree_p = drop.tip(tree, setdiff(tree$tip.label, align_p$matings.especies))
#nicho_p = nicho_f[35:53,]

#
df = data.frame(mds_p$especies, mds_p$dados.PARVORDER, mds_p$dados.SOCIAL_ORGANIZATION,
                mds_p$dados.MATING_SYSTEM, mds_p$dados.DOMINANCE,
                mds_p$dados.AGGRESSION, evolvas_p$Dimorfism, 
                evolvas_p$Integration, evolvas_p$Size, mds_p$dados.PROP_MALES_FEMALES,
                mds_p$scores.fit....1., mds_p$scores.fit....2., align_p$normas)
                #nicho_p$climacv, nicho_p$CVn)


df <- df[match(tree_p$tip.label, df$mds_p.especies), ]

# Modelos
# BM
df_bm <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas
)

df_bm <- df_bm[match(tree_p$tip.label, df$mds_p.especies), ]

sygma = seq(0.002, 0.007, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_p,
  response = df_bm$response,
  species = df_bm$species,
  sigma2_y_values = sygma
)
plot(fit_bm)

#save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.03, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_p,
  response = df_bm$response,
  species =  df_bm$species,
  a_values = alpha_vals
)
plot(fit_ou)

#save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU.RData")

# OU Size
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas,
  size = scale(df$evolvas_p.Size)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.especies), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_size <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$size,
  a_values = alpha_vals
)
plot(fit_ou_size)

#save(fit_ou_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_size.RData")

# OU NMDS1
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.scores.fit....1.)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.especies), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds1 <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds1,
  a_values = alpha_vals
)
plot(fit_ou_nmds1)

#save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds2 = scale(df$mds_p.scores.fit....2.)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.especies), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds2 <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds2
)
plot(fit_ou_nmds2)

#save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS2.RData")

# OU Integration
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas,  
  integration = scale(df$evolvas_p.Integration)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.especies), ]

alpha_vals <- seq(0.01, 0.07, length.out = 20)
fit_ou_integration <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("integration")],
  a_values = alpha_vals
)
plot(fit_ou_integration)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Integration.RData")

# # OU Clima
# df_ou <- data.frame(
#   species = factor(df$mds_p.especies, levels = tree_p$tip.label),
#   response = align_p$normas,
#   clima = df$nicho_p.climacv
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$species), ]
# 
# fit_ou_clima <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou$clima
# )
# summary(fit_ou_clima)
# 
# #save(fit_ou_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Clima.RData")
# 
# # OU CVn
# df_ou <- data.frame(
#   species = factor(df$mds_p.especies, levels = tree_p$tip.label),
#   response = align_p$normas,
#   cvn = df$nicho_p.CVn
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$species), ]
# 
# alpha_vals <- seq(0.01, 0.1, length.out = 20)
# fit_ou_cvn <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou$cvn,
#   a_values = alpha_vals
# )
# plot(fit_ou_cvn)
# 
# #save(fit_ou_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_CVn.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.scores.fit....1.),
  nmds2 = scale(df$mds_p.scores.fit....2.)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.especies), ]

fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2")],
)
summary(fit_ou_nmds1_nmds2)

#save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Integration
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.scores.fit....1.),
  nmds2 = scale(df$mds_p.scores.fit....2.),
  integration = scale(df$evolvas_p.Integration)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.especies), ]

fit_ou_nmds1_nmds2_integration <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "integration")]
)
summary(fit_ou_nmds1_nmds2_integration)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU NMDS1 + NMDS2 + Integration + Size
df_ou <- data.frame(
  species = factor(df$mds_p.especies, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.scores.fit....1.),
  nmds2 = scale(df$mds_p.scores.fit....2.),
  integration = scale(df$evolvas_p.Integration),
  size = scale(df$evolvas_p.Size)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.especies), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds1_nmds2_integration_size <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size")],
  a_values = alpha_vals
)
plot(fit_ou_nmds1_nmds2_integration_size)

#save(fit_ou_nmds1_nmds2_integration_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size.RData")

# # OU NMDS1 + NMDS2 + Integration + Size + Clima
# df_ou <- data.frame(
#   species = factor(df$mds_p.especies, levels = tree_p$tip.label),
#   response = align_p$normas,
#   nmds1 = df$mds_p.scores.fit....1.,
#   nmds2 = df$mds_p.scores.fit....2.,
#   integration = df$evolvas_p.Integration,
#   size = df$evolvas_p.Size,
#   clima = df$nicho_p.climacv
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$species), ]
# 
# alpha_vals <- seq(0.01, 0.1, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima)
# 
# #save(fit_ou_nmds1_nmds2_integration_size_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima.RData")
# 
# # OU NMDS1 + NMDS2 + Integration + Size + Clima + CVn
# df_ou <- data.frame(
#   species = factor(df$mds_p.especies, levels = tree_p$tip.label),
#   response = align_p$normas,
#   nmds1 = df$mds_p.scores.fit....1.,
#   nmds2 = df$mds_p.scores.fit....2.,
#   integration = df$evolvas_p.Integration,
#   size = df$evolvas_p.Size,
#   clima = df$nicho_p.climacv,
#   cvn = df$nicho_p.CVn
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$species), ]
# 
# alpha_vals <- seq(0.01, 0.1, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima_cvn <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima", "cvn")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima_cvn)
# 
# #save(fit_ou_nmds1_nmds2_integration_size_clima_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima_CVn.RData")

fits = list(
  BM = fit_bm,
  OU = fit_ou,
  OU_SIZE = fit_ou_size,
  OU_NMDS1 = fit_ou_nmds1,
  OU_NMDS2 = fit_ou_nmds2,
  OU_INTEGRATION = fit_ou_integration,
  #OU_CLIMA = fit_ou_clima,
  #OU_CVN = fit_ou_cvn,
  OU_NMDS1_NMDS2 = fit_ou_nmds1_nmds2,
  OU_NMDS1_NMDS2_INTEGRATION = fit_ou_nmds1_nmds2_integration,
  OU_NMDS1_NMDS2_INTEGRATION_SIZE = fit_ou_nmds1_nmds2_integration_size
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fit_ou_nmds1_nmds2_integration_size_clima,
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fit_ou_nmds1_nmds2_integration_size_clima_cvn
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fits.RData")

#################################################################################
# Haplorrhini
#
df = data.frame(mds_f$especies, mds_f$dados.PARVORDER, mds_f$dados.SOCIAL_ORGANIZATION,
                mds_f$dados.MATING_SYSTEM, mds_f$dados.AGGRESSION, evolvas_f$Dimorfism, 
                evolvas_f$Integration, evolvas_f$Size, mds_f$dados.PROP_MALES_FEMALES,
                mds_f$scores.fit....1., mds_f$scores.fit....2., align_f$normas)
                #, nicho_f$climacv, nicho_f$CVn)


df <- df[match(tree$tip.label, df$mds_f.especies), ]

# Modelos
# BM
df_bm <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas
)

df_bm <- df_bm[match(tree$tip.label, df$mds_f.especies), ]

sygma = seq(0.002, 0.007, length.out = 20)
fit_bm <- brown.fit(
  phy = tree,
  response = df_bm$response,
  species = df_bm$species,
  sigma2_y_values = sygma
)
plot(fit_bm)

#save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.03, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree,
  response = df_bm$response,
  species =  df_bm$species,
  a_values = alpha_vals
)
plot(fit_ou)

#save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU.RData")

# OU Size
df_ou <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas,
  size = scale(df$evolvas_f.Size)
)

df_ou <- df_ou[match(tree$tip.label, df$mds_f.especies), ]

fit_ou_size <- slouch.fit(
  phy = tree,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$size
)
summary(fit_ou_size)

#save(fit_ou_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_size.RData")

# OU NMDS1
df_ou <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas,
  nmds1 = scale(df$mds_f.scores.fit....1.)
)

df_ou <- df_ou[match(tree$tip.label, df$mds_f.especies), ]

fit_ou_nmds1 <- slouch.fit(
  phy = tree,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds1
)
summary(fit_ou_nmds1)

#save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas,
  nmds2 = scale(df$mds_f.scores.fit....2.)
)

df_ou <- df_ou[match(tree$tip.label, df$mds_f.especies), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds2 <- slouch.fit(
  phy = tree,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou$nmds2
)
plot(fit_ou_nmds2)

#save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS2.RData")

# OU Integration
df_ou <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas,  
  integration = scale(df$evolvas_f.Integration)
)

df_ou <- df_ou[match(tree$tip.label, df$mds_f.especies), ]

fit_ou_integration <- slouch.fit(
  phy = tree,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("integration")]
)
summary(fit_ou_integration)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Integration.RData")

# # OU Clima 
# df_ou <- data.frame(
#   species = factor(df$mds_f.especies, levels = tree$tip.label),
#   response = align_f$normas,
#   clima = df$nicho_f.climacv
# )
# 
# df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]
# 
# fit_ou_clima <- slouch.fit(
#   phy = tree,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou$clima
#   )
# summary(fit_ou_clima)
# 
# #save(fit_ou_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Clima.RData")
# 
# # OU CVn
# df_ou <- data.frame(
#   species = factor(df$mds_f.especies, levels = tree$tip.label),
#   response = align_f$normas,
#   cvn = df$nicho_f.CVn
# )
# 
# df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]
# 
# fit_ou_cvn <- slouch.fit(
#   phy = tree,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou$cvn
# )
# summary(fit_ou_cvn)
# 
# #save(fit_ou_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_CVn.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas,
  nmds1 = scale(df$mds_f.scores.fit....1.),
  nmds2 = scale(df$mds_f.scores.fit....2.)
)

df_ou <- df_ou[match(tree$tip.label, df$mds_f.especies), ]

fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2")]
)
plot(fit_ou_nmds1_nmds2)

#save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Integration
df_ou <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas,
  nmds1 = scale(df$mds_f.scores.fit....1.),
  nmds2 = scale(df$mds_f.scores.fit....2.),
  integration = scale(df$evolvas_f.Integration)
)

df_ou <- df_ou[match(tree$tip.label, df$mds_f.especies), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds1_nmds2_integration <- slouch.fit(
  phy = tree,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "integration")],
  a_values = alpha_vals
)
plot(fit_ou_nmds1_nmds2_integration)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU NMDS1 + NMDS2 + Integration + Size
df_ou <- data.frame(
  species = factor(df$mds_f.especies, levels = tree$tip.label),
  response = df$align_f.normas,
  nmds1 = scale(df$mds_f.scores.fit....1.),
  nmds2 = scale(df$mds_f.scores.fit....2.),
  integration = scale(df$evolvas_f.Integration),
  size = scale(df$evolvas_f.Size)
)

df_ou <- df_ou[match(tree$tip.label, df$mds_f.especies), ]

fit_ou_nmds1_nmds2_integration_size <- slouch.fit(
  phy = tree,
  response = df_ou$response,
  species =  df_ou$species,
  random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size")]
)
summary(fit_ou_nmds1_nmds2_integration_size)

#save(fit_ou_nmds1_nmds2_integration_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size.RData")

# # OU NMDS1 + NMDS2 + Integration + Size + Clima
# df_ou <- data.frame(
#   species = factor(df$mds_f.especies, levels = tree$tip.label),
#   response = align_f$normas,
#   nmds1 = df$mds_f.scores.fit....1.,
#   nmds2 = df$mds_f.scores.fit....2.,
#   integration = df$evolvas_f.Integration,
#   size = df$evolvas_f.Size,
#   clima = df$nicho_f.climacv
# )
# 
# df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]
# 
# fit_ou_nmds1_nmds2_integration_size_clima <- slouch.fit(
#   phy = tree,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima")]
# )
# summary(fit_ou_nmds1_nmds2_integration_size_clima)
# 
# #save(fit_ou_nmds1_nmds2_integration_size_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima.RData")
# 
# # OU NMDS1 + NMDS2 + Integration + Size + Clima + CVn
# df_ou <- data.frame(
#   species = factor(df$mds_f.especies, levels = tree$tip.label),
#   response = align_f$normas,
#   nmds1 = df$mds_f.scores.fit....1.,
#   nmds2 = df$mds_f.scores.fit....2.,
#   integration = df$evolvas_f.Integration,
#   size = df$evolvas_f.Size,
#   clima = df$nicho_f.climacv,
#   cvn = df$nicho_f.CVn
# )
# 
# df_ou <- df_ou[match(tree$tip.label, df_ou$species), ]
# 
# fit_ou_nmds1_nmds2_integration_size_clima_cvn <- slouch.fit(
#   phy = tree,
#   response = df_ou$response,
#   species =  df_ou$species,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima", "cvn")]
# )
# summary(fit_ou_nmds1_nmds2_integration_size_clima_cvn)
# 
# save(fit_ou_nmds1_nmds2_integration_size_clima_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/haplorrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima_CVn.RData")

fits = list(
  BM = fit_bm,
  OU = fit_ou,
  OU_SIZE = fit_ou_size,
  OU_NMDS1 = fit_ou_nmds1,
  OU_NMDS2 = fit_ou_nmds2,
  OU_INTEGRATION = fit_ou_integration,
  #OU_CLIMA = fit_ou_clima,
  #OU_CVN = fit_ou_cvn,
  OU_NMDS1_NMDS2 = fit_ou_nmds1_nmds2,
  OU_NMDS1_NMDS2_INTEGRATION = fit_ou_nmds1_nmds2_integration,
  OU_NMDS1_NMDS2_INTEGRATION_SIZE = fit_ou_nmds1_nmds2_integration_size
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA = fit_ou_nmds1_nmds2_integration_size_clima,
  #OU_NMDS1_NMDS2_INTEGRATION_SIZE_CLIMA_CVN = fit_ou_nmds1_nmds2_integration_size_clima_cvn
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/haplorrhini_fits.RData")
