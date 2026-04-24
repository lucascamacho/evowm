# SLOUCH POR GENERO
# CATARRHINI, PLATYRRHINI E HAPLORRHYNI

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

library(slouch)
library(ape)
library(dplyr)

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")
#nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
#  nicho$genus,
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
plot(tree2)
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

################################################################################
# A PARTIR DAQUI
# Catarrhini 1:16 com nicho, 1:21 sem nicho
align_c = align_f[1:21,]
evolvas_c = evolvas_f[1:21,]
mds_c = mds_f[1:21,]
tree_c = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_c$matings.genus))
#nicho_c = nicho_f[1:16,]

#
df = data.frame(mds_c$genus, mds_c$dados.PARVORDER, mds_c$dados.SOCIAL_ORGANIZATION,
                mds_c$dados.AGGRESSION, mds_c$dados.DOMINANCE, mds_c$dados.ALL_MALE_GROUPS,
                mds_c$dados.FURTIVE_COPULATION, evolvas_c$Dimorfism, 
                evolvas_c$Integration, evolvas_c$Size, mds_c$vegan..scores.fit..sites...1., 
                mds_c$vegan..scores.fit..sites...2., align_c$normas)
#                , nicho_c$climacv_mean, nicho_c$CVn_mean)

df <- df[match(tree_c$tip.label, df$mds_c.genus), ]

# Modelos
# BM
df_bm <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = align_c$normas
)
df_bm <- df_bm[match(tree_c$tip.label, df_bm$genus), ]

sygma = seq(0.01, 0.05, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_c,
  response = df_bm$response,
  species = df_bm$genus,
  sigma2_y_values = sygma
)
plot(fit_bm)

#save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.3, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_c,
  response = df_bm$response,
  species =  df_bm$genus,
  a_values = alpha_vals
)
plot(fit_ou)

#save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU.RData")

# OU Size
df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = df$align_c.normas,
  size = scale(df$evolvas_c.Size)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

fit_ou_size <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$size
)
plot(fit_ou_size)

#save(fit_ou_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_size.RData")

# OU NMDS1
df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.vegan..scores.fit..sites...1.)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

alpha_vals <- seq(0.01, 0.4, length.out = 20)
fit_ou_nmds1 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$nmds1,
  a_values = alpha_vals
)
plot(fit_ou_nmds1)

#save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds2 = scale(df$mds_c.vegan..scores.fit..sites...2.)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

alpha_vals <- seq(0.01, 0.2, length.out = 20)
fit_ou_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$nmds2,
  a_values = alpha_vals
)
plot(fit_ou_nmds2)

#save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS2.RData")

# OU Integration 
df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = align_c$normas,
  integration = scale(df$evolvas_c.Integration)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$integration,
  a_values = alpha_vals
)
plot(fit_ou_integration)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Integration.RData")

# OU Clima
# df_ou <- data.frame(
#   genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
#   response = align_c$normas,
#   clima = df$nicho_c.climacv_mean
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.9, length.out = 20)
# fit_ou_clima <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou$clima,
#   a_values = alpha_vals
# )
# plot(fit_ou_clima)

#save(fit_ou_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Clima.RData")

# OU CVn
# df_ou <- data.frame(
#   genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
#   response = align_c$normas,
#   cvn = df$nicho_c.CVn_mean
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(100, 10000, length.out = 20)
# fit_ou_cvn <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou$cvn,
#   a_values = alpha_vals
# )
# summary(fit_ou_cvn)

#save(fit_ou_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/catarrhini_fit_OU_CVn.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_c.vegan..scores.fit..sites...2.)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

alpha_vals <- seq(0.01, 0.2, length.out = 20)
fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2")],
  a_values = alpha_vals
)
plot(fit_ou_nmds1_nmds2)

#save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Integration
df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_c.vegan..scores.fit..sites...2.),
  integration = scale(df$evolvas_c.Integration)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

fit_ou_nmds1_nmds2_integration <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2", "integration")]
)
summary(fit_ou_nmds1_nmds2_integration)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU NMDS1 + NMDS2 + Integration + Size
df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  response = align_c$normas,
  nmds1 = scale(df$mds_c.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_c.vegan..scores.fit..sites...2.),
  integration = scale(df$evolvas_c.Integration),
  size = scale(df$evolvas_c.Size)
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

fit_ou_nmds1_nmds2_integration_size <- slouch.fit(
  phy = tree_c,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2", "integration", "size")]
)
summary(fit_ou_nmds1_nmds2_integration_size)

# OU NMDS1 + NMDS2 + Integration + Size + Clima
# df_ou <- data.frame(
#   genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
#   response = align_c$normas,
#   nmds1 = df$mds_c.vegan..scores.fit..sites...1.,
#   nmds2 = df$mds_c.vegan..scores.fit..sites...2.,
#   integration = df$evolvas_c.Integration,
#   size = df$evolvas_c.Size,
#   clima = df$nicho_c.climacv_mean
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.9, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima)

#save(fit_ou_nmds1_nmds2_integration_size_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima.RData")

# OU NMDS1 + NMDS2 + Integration + Size + Clima + CVn
# df_ou <- data.frame(
#   genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
#   response = align_c$normas,
#   nmds1 = df$mds_c.vegan..scores.fit..sites...1.,
#   nmds2 = df$mds_c.vegan..scores.fit..sites...2.,
#   integration = df$evolvas_c.Integration,
#   size = df$evolvas_c.Size,
#   clima = df$nicho_c.climacv_mean,
#   cvn = df$nicho_c.CVn_mean
# )
# 
# df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.9, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima_cvn <- slouch.fit(
#   phy = tree_c,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima", "cvn")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima_cvn)

#save(fit_ou_nmds1_nmds2_integration_size_clima_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima_CVn.RData")

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

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/catarrhini_fits.RData")

#################################################################################
# Platyrrhini 22:37 sem nicho, 17:28 com nicho
align_p = align_f[22:37,]
evolvas_p = evolvas_f[22:37,]
mds_p = mds_f[22:37,]
tree_p = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_p$matings.genus))
#nicho_p = nicho_f[17:28,]

#
df = data.frame(mds_p$genus, mds_p$dados.PARVORDER, mds_p$dados.SOCIAL_ORGANIZATION,
                mds_p$dados.MATING_SYSTEM, mds_p$dados.AGGRESSION, evolvas_p$Dimorfism, 
                evolvas_p$Integration, evolvas_p$Size, mds_p$dados.PROP_MALES_FEMALES,
                mds_p$vegan..scores.fit..sites...1., mds_p$vegan..scores.fit..sites...2., 
                align_p$normas) 
#                , nicho_p$climacv_mean, nicho_p$CVn_mean)


df <- df[match(tree_p$tip.label, df$mds_p.genus), ]

# Modelos
# BM
df_bm <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas
)

df_bm <- df_bm[match(tree_p$tip.label, df$mds_p.genus), ]

sygma = seq(0.002, 0.007, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_p,
  response = df_bm$response,
  species = df_bm$genus,
  sigma2_y_values = sygma
)
plot(fit_bm)

#save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.3, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_p,
  response = df_bm$response,
  species =  df_bm$genus,
  a_values = alpha_vals
)
plot(fit_ou)

#save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU.RData")

# OU Size 
df_ou <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas,
  size = scale(df$evolvas_p.Size)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.genus), ]

fit_ou_size <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$size,
)
plot(fit_ou_size)

#save(fit_ou_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_size.RData")

# OU NMDS1
df_ou <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.vegan..scores.fit..sites...1.)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.genus), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds1 <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$nmds1,
  a_values = alpha_vals
)
plot(fit_ou_nmds1)

#save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds2 = scale(df$mds_p.vegan..scores.fit..sites...2.)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.genus), ]

alpha_vals <- seq(0.01, 0.1, length.out = 20)
fit_ou_nmds2 <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$nmds2
)
plot(fit_ou_nmds2)

#save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS2.RData")

# OU Integration
df_ou <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas,  
  integration = scale(df$evolvas_p.Integration)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.genus), ]

alpha_vals <- seq(10, 100, length.out = 20)
fit_ou_integration <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("integration")],
  a_values = alpha_vals
)
summary(fit_ou_integration)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/platyrrhini_fit_OU_Integration.RData")

# OU Clima
# df_ou <- data.frame(
#   genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
#   response = align_p$normas,
#   clima = df$nicho_p.climacv_mean
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.1, length.out = 20)
# fit_ou_clima <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou$clima,
#   a_values = alpha_vals
# )
# plot(fit_ou_clima)

#save(fit_ou_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_Clima.RData")

# OU CVn
# df_ou <- data.frame(
#   genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
#   response = align_p$normas,
#   cvn = df$nicho_p.CVn_mean
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.1, length.out = 20)
# fit_ou_cvn <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou$cvn,
#   a_values = alpha_vals
# )
# plot(fit_ou_cvn)

#save(fit_ou_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_CVn.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_p.vegan..scores.fit..sites...2.)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.genus), ]

fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2")],
)
plot(fit_ou_nmds1_nmds2)

#save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Integration
df_ou <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_p.vegan..scores.fit..sites...2.),
  integration = scale(df$evolvas_p.Integration)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.genus), ]

fit_ou_nmds1_nmds2_integration <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2", "integration")]
)
summary(fit_ou_nmds1_nmds2_integration)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU NMDS1 + NMDS2 + Integration + Size
df_ou <- data.frame(
  genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
  response = df$align_p.normas,
  nmds1 = scale(df$mds_p.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_p.vegan..scores.fit..sites...2.),
  integration = scale(df$evolvas_p.Integration),
  size = scale(df$evolvas_p.Size)
)

df_ou <- df_ou[match(tree_p$tip.label, df$mds_p.genus), ]

alpha_vals <- seq(0.1, 0.9, length.out = 20)
fit_ou_nmds1_nmds2_integration_size <- slouch.fit(
  phy = tree_p,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2", "integration", "size")],
  a_values = alpha_vals
)
plot(fit_ou_nmds1_nmds2_integration_size)

#save(fit_ou_nmds1_nmds2_integration_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size.RData")

# OU NMDS1 + NMDS2 + Integration + Size + Clima
# df_ou <- data.frame(
#   genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
#   response = align_p$normas,
#   nmds1 = df$mds_p.vegan..scores.fit..sites...1.,
#   nmds2 = df$mds_p.vegan..scores.fit..sites...2.,
#   integration = df$evolvas_p.Integration,
#   size = df$evolvas_p.Size,
#   clima = df$nicho_p.climacv_mean
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.5, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima)

#save(fit_ou_nmds1_nmds2_integration_size_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima.RData")

# OU NMDS1 + NMDS2 + Integration + Size + Clima + CVn
# df_ou <- data.frame(
#   genus = factor(df$mds_p.genus, levels = tree_p$tip.label),
#   response = align_p$normas,
#   nmds1 = df$mds_p.vegan..scores.fit..sites...1.,
#   nmds2 = df$mds_p.vegan..scores.fit..sites...2.,
#   integration = df$evolvas_p.Integration,
#   size = df$evolvas_p.Size,
#   clima = df$nicho_p.climacv_mean,
#   cvn = df$nicho_p.CVn_mean
# )
# 
# df_ou <- df_ou[match(tree_p$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.9, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima_cvn <- slouch.fit(
#   phy = tree_p,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima", "cvn")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima_cvn)

#save(fit_ou_nmds1_nmds2_integration_size_clima_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima_CVn.RData")

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

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/platyrrhini_fits.RData")

#################################################################################
# Haplorrhini
df = data.frame(mds_f$genus, mds_f$dados.PARVORDER, mds_f$dados.SOCIAL_ORGANIZATION,
                mds_f$dados.MATING_SYSTEM, mds_f$dados.AGGRESSION, evolvas_f$Dimorfism, 
                evolvas_f$Integration, 
                evolvas_f$Size, 
                mds_f$vegan..scores.fit..sites...1., 
                mds_f$vegan..scores.fit..sites...2., align_f$normas)
#                , nicho_f$climacv_mean, nicho_f$CVn_mean)

df <- df[match(tree_genus$tip.label, df$mds_f.genus), ]

# Modelos
# BM
df_bm <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = align_f$normas
)
df_bm <- df_bm[match(tree_genus$tip.label, df_bm$genus), ]

sygma = seq(0.01, 0.05, length.out = 20)
fit_bm <- brown.fit(
  phy = tree_genus,
  response = df_bm$response,
  species = df_bm$genus,
  sigma2_y_values = sygma
)
plot(fit_bm)

#save(fit_bm, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_BM.RData")

# OU SEM PREDITORES
alpha_vals <- seq(0.01, 0.3, length.out = 20)
fit_ou <- slouch.fit(
  phy = tree_genus,
  response = df_bm$response,
  species =  df_bm$genus,
  a_values = alpha_vals
)
plot(fit_ou)

#save(fit_ou, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU.RData")

# OU Size
df_ou <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = df$align_f.normas,
  size = scale(df$evolvas_f.Size)
)

df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]

fit_ou_size <- slouch.fit(
  phy = tree_genus,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$size
)
summary(fit_ou_size)

#save(fit_ou_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_size.RData")

# OU NMDS1
df_ou <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = align_f$normas,
  nmds1 = scale(df$mds_f.vegan..scores.fit..sites...1.)
)

df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]

alpha_vals <- seq(0.01, 0.4, length.out = 20)
fit_ou_nmds1 <- slouch.fit(
  phy = tree_genus,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$nmds1,
  a_values = alpha_vals
)
plot(fit_ou_nmds1)

#save(fit_ou_nmds1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1.RData")

# OU NMDS2
df_ou <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = align_f$normas,
  nmds2 = scale(df$mds_f.vegan..scores.fit..sites...2.)
)

df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]

fit_ou_nmds2 <- slouch.fit(
  phy = tree_genus,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$nmds2
)
plot(fit_ou_nmds2)

#save(fit_ou_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS2.RData")

# OU Integration 
df_ou <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = align_f$normas,
  integration = scale(df$evolvas_f.Integration)
)

df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]

fit_ou_integration <- slouch.fit(
  phy = tree_genus,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou$integration
)

summary(fit_ou_integration)

#save(fit_ou_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Integration.RData")

# # OU Clima
# df_ou <- data.frame(
#   genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
#   response = align_f$normas,
#   clima = df$nicho_f.climacv_mean
# )
# 
# df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.9, length.out = 20)
# fit_ou_clima <- slouch.fit(
#   phy = tree_genus,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou$clima,
#   a_values = alpha_vals
# )
# plot(fit_ou_clima)
# 
# #save(fit_ou_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_Clima.RData")
# 
# # OU CVn
# df_ou <- data.frame(
#   genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
#   response = align_f$normas,
#   cvn = df$nicho_f.CVn_mean
# )
# 
# df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]
# 
# fit_ou_cvn <- slouch.fit(
#   phy = tree_genus,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou$cvn
# )
# summary(fit_ou_cvn)

#save(fit_ou_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_CVn.RData")

# OU NMDS1 + NMDS2
df_ou <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = align_f$normas,
  nmds1 = scale(df$mds_f.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_f.vegan..scores.fit..sites...2.)
)

df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]

alpha_vals <- seq(0.01, 0.2, length.out = 20)
fit_ou_nmds1_nmds2 <- slouch.fit(
  phy = tree_genus,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2")],
  a_values = alpha_vals
)
plot(fit_ou_nmds1_nmds2)

#save(fit_ou_nmds1_nmds2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2.RData")

# OU NMDS1 + NMDS2 + Integration
df_ou <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = align_f$normas,
  nmds1 = scale(df$mds_f.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_f.vegan..scores.fit..sites...2.),
  integration = scale(df$evolvas_f.Integration)
)

df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]

fit_ou_nmds1_nmds2_integration <- slouch.fit(
  phy = tree_genus,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2", "integration")]
)
summary(fit_ou_nmds1_nmds2_integration)

#save(fit_ou_nmds1_nmds2_integration, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration.RData")

# OU NMDS1 + NMDS2 + Integration + Size
df_ou <- data.frame(
  genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
  response = align_f$normas,
  nmds1 = scale(df$mds_f.vegan..scores.fit..sites...1.),
  nmds2 = scale(df$mds_f.vegan..scores.fit..sites...2.),
  integration = scale(df$evolvas_f.Integration),
  size = scale(df$evolvas_f.Size)
)

df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]

fit_ou_nmds1_nmds2_integration_size <- slouch.fit(
  phy = tree_genus,
  response = df_ou$response,
  species =  df_ou$genus,
  direct.cov = df_ou[, c("nmds1", "nmds2", "integration", "size")]
)
summary(fit_ou_nmds1_nmds2_integration_size)

#save(fit_ou_nmds1_nmds2_integration_size, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/haplorrhini_fit_OU_NMDS1_NMDS2_Integration_Size.RData")

# # OU NMDS1 + NMDS2 + Integration + Size + Clima
# df_ou <- data.frame(
#   genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
#   response = align_f$normas,
#   nmds1 = df$mds_f.vegan..scores.fit..sites...1.,
#   nmds2 = df$mds_f.vegan..scores.fit..sites...2.,
#   integration = df$evolvas_f.Integration,
#   size = df$evolvas_f.Size,
#   clima = df$nicho_f.climacv_mean
# )
# 
# df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.9, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima <- slouch.fit(
#   phy = tree_genus,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima)
# 
# #save(fit_ou_nmds1_nmds2_integration_size_clima, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima.RData")
# 
# # OU NMDS1 + NMDS2 + Integration + Size + Clima + CVn
# df_ou <- data.frame(
#   genus = factor(df$mds_f.genus, levels = tree_genus$tip.label),
#   response = align_f$normas,
#   nmds1 = df$mds_f.vegan..scores.fit..sites...1.,
#   nmds2 = df$mds_f.vegan..scores.fit..sites...2.,
#   integration = df$evolvas_f.Integration,
#   size = df$evolvas_f.Size,
#   clima = df$nicho_f.climacv_mean,
#   cvn = df$nicho_f.CVn_mean
# )
# 
# df_ou <- df_ou[match(tree_genus$tip.label, df_ou$genus), ]
# 
# alpha_vals <- seq(0.01, 0.9, length.out = 20)
# fit_ou_nmds1_nmds2_integration_size_clima_cvn <- slouch.fit(
#   phy = tree_genus,
#   response = df_ou$response,
#   species =  df_ou$genus,
#   random.cov = df_ou[, c("nmds1", "nmds2", "integration", "size", "clima", "cvn")],
#   a_values = alpha_vals
# )
# plot(fit_ou_nmds1_nmds2_integration_size_clima_cvn)
# 
# #save(fit_ou_nmds1_nmds2_integration_size_clima_cvn, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/haplorrhini_fit_OU_NMDS1_NMDS2_Integration_Size_Clima_CVn.RData")

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

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/haplorrhini_fits.RData")
