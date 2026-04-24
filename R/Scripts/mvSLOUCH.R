setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(mvSLOUCH)
library(ape)
library(PCMBaseCpp)

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
tree = ape::read.nexus(filename)
species = mds$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

align_filtrado <- align[align$matings.especies %in% mds$especies, ]
evolvas_filtrado <- evolvas[evolvas$species %in% mds$especies, ]

species62 = align_filtrado$matings.especies

by_trait <- medias$ByTrait_Averages

by_trait <- by_trait[names(by_trait) %in% species62]

# inicializar lista ou matriz pra guardar o dimorfismo
dimorfismo_list <- list()

# loop sobre espécies
for(sp in names(by_trait)) {
  macho <- by_trait[[sp]]$Machos       # vetor com 39 traços
  femea <- by_trait[[sp]]$Fêmeas       # vetor com 39 traços
  
  # calcular dimorfismo normalizado
  dimorfismo_list[[sp]] <- (macho - femea) / ((geomean(macho) + geomean(femea)) / 2)
}

dimorfismo_mat <- do.call(rbind, dimorfismo_list)
rownames(dimorfismo_mat) <- names(dimorfismo_list)
dimorfismo_mat <- dimorfismo_mat[tree$tip.label, ]

# VARIÁVEIS ORIGINAIS
# Catarrhini: SOCIAL_ORG, MATING SYSTEM(quase), AGGRESSION

#####################################################
# A PARTIR DAQUI
# Catarrhini
align_filtrado_c = align_filtrado[1:38,]
evolvas_filtrado_c = evolvas_filtrado[1:38,]
mds_c = mds[1:38,]
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_filtrado_c$matings.especies))
dimorfismo_mat_c = dimorfismo_mat[1:38,]

#
df = data.frame(mds_c$especies, mds_c$dados.PARVORDER, mds_c$dados.SOCIAL_ORGANIZATION,
                mds_c$dados.MATING_SYSTEM, mds_c$dados.PROP_MALES_FEMALES, mds_c$dados.DOMINANCE,
                mds_c$dados.AGGRESSION, mds_c$dados.ALL_MALE_GROUPS, mds_c$dados.FURTIVE_COPULATION,
                mds_c$dados.INFANTICIDE, mds_c$dados.MULTILEVEL_SOCIETY, mds_c$dados.N_PAPERS,
                align_filtrado_c$dimor, align_filtrado_c$normas, evolvas_filtrado_c$Size, evolvas_filtrado_c$Integration)

df <- df[match(tree_c$tip.label, df$mds_c.especies), ]

# Modelos
# BM
fit_BM <- BrownianMotionModel(
  phyltree = tree_c,
  mData = dimorfismo_mat_c
)
save(fit_BM, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_BM.RData")

# OU Simple
fit_OU <- ouchModel(
  phyltree = tree_c,
  mData = dimorfismo_mat_c,
  Atype = "Diagonal",        
  estimate.root.state = FALSE
)

save(fit_OU, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_SIMPLE_OU.RData")

# OU SOCIAL ORG
mData_OU_Social <- cbind(dimorfismo_mat_c, df$mds_c.dados.SOCIAL_ORGANIZATION)
colnames(mData_OU_Social)[40] <- "SOCIAL"
rownames(mData_OU_Social) <- df$mds_c.especies
mData_OU_Social <- mData_OU_Social[tree_c$tip.label, ]

fit_OU_SOCIAL <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_Social,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_SOCIAL.RData")

# OU MATING SYSTEM
mData_OU_Mating <- cbind(dimorfismo_mat_c, df$mds_c.dados.MATING_SYSTEM)
colnames(mData_OU_Mating)[40] <- "MATING"
rownames(mData_OU_Mating) <- df$mds_c.especies
mData_OU_Mating <- mData_OU_Mating[tree_c$tip.label, ]

fit_OU_MATINGS <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_Mating,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_MATINGS, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_MATINGS.RData")

# OU AGGRESSION
mData_OU_Agression <- cbind(dimorfismo_mat_c, df$mds_c.dados.AGGRESSION)
colnames(mData_OU_Agression)[40] <- "AGGRESSION"
rownames(mData_OU_Agression) <- df$mds_c.especies
mData_OU_Agression <- mData_OU_Agression[tree_c$tip.label, ]

fit_OU_AGGRESSION <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_Agression,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_AGGRESSION, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_AGGRESSION.RData")

# OU INTEGRATION
mData_OU_Integration <- cbind(dimorfismo_mat_c, df$evolvas_filtrado_c.Integration)
colnames(mData_OU_Integration)[40] <- "INTEGRATION"
rownames(mData_OU_Integration) <- df$mds_c.especies
mData_OU_Integration <- mData_OU_Integration[tree_c$tip.label, ]

fit_OU_INTEGRATION <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_Integration,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_INTEGRATION, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_INTEGRATION.RData")

# OU SOCIAL + MATINGS
mData_OU_Social_Matings <- cbind(dimorfismo_mat_c, 
                              df$mds_c.dados.SOCIAL_ORGANIZATION,
                              df$mds_c.dados.MATING_SYSTEM)
colnames(mData_OU_Social_Matings)[40:41] <- c("SOCIAL", "MATINGS")
rownames(mData_OU_Social_Matings) <- df$mds_c.especies
mData_OU_Social_Matings <- mData_OU_Social_Matings[tree_c$tip.label, ]

fit_OU_SOCIAL_MATINGS <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_Social_Matings,
  kY = 39,
  predictors = c(40, 41),
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL_MATINGS, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_SOCIAL_MATINGS.RData")

# OU SOCIAL + MATINGS + AGGRESSION
mData_OU_Social_Matings_Aggression <- cbind(dimorfismo_mat_c,
                                 df$mds_c.dados.SOCIAL_ORGANIZATION,
                                 df$mds_c.dados.MATING_SYSTEM,
                                 df$mds_c.dados.AGGRESSION)
colnames(mData_OU_Social_Matings_Aggression)[40:42] <- c("SOCIAL", "MATINGS", "AGGRESSION")
rownames(mData_OU_Social_Matings_Aggression) <- df$mds_c.especies
mData_OU_Social_Matings_Aggression <- mData_OU_Social_Matings_Aggression[tree_c$tip.label, ]

fit_OU_SOCIAL_MATINGS_AGGRESSION <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_Social_Matings_Aggression,
  kY = 39,
  predictors = c(40, 42),
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL_MATINGS_AGGRESSION, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_SOCIAL_MATINGS_AGGRESSION.RData")

# OU SOCIAL + MATINGS + AGGRESSION + INTEGRATION
mData_OU_Social_Matings_Aggression_Integration <- cbind(dimorfismo_mat_c,
                                            df$mds_c.dados.SOCIAL_ORGANIZATION,
                                            df$mds_c.dados.MATING_SYSTEM,
                                            df$mds_c.dados.AGGRESSION,
                                            df$evolvas_filtrado_c.Integration)
colnames(mData_OU_Social_Matings_Aggression_Integration)[40:43] <- c("SOCIAL", "MATINGS", "AGGRESSION", "INTEGRATION")
rownames(mData_OU_Social_Matings_Aggression_Integration) <- df$mds_c.especies
mData_OU_Social_Matings_Aggression_Integration <- mData_OU_Social_Matings_Aggression_Integration[tree_c$tip.label, ]

fit_OU_SOCIAL_MATINGS_AGGRESSION_INTEGRATION <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_Social_Matings_Aggression_Integration,
  kY = 39,
  predictors = c(40, 43),
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL_MATINGS_AGGRESSION_INTEGRATION, 
     file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_SOCIAL_MATINGS_AGGRESSION_INTEGRATION.RData")

fits = list(
  BM = fit_BM,
  OU_SIMPLE = fit_OU,
  OU_SOCIAL = fit_OU_SOCIAL,
  OU_MATINGS = fit_OU_MATINGS,
  OU_AGGRESSION = fit_OU_AGGRESSION,
  OU_INTEGRATION = fit_OU_INTEGRATION,
  OU_SOCIAL_MATINGS = fit_OU_SOCIAL_MATINGS,
  OU_SOCIAL_MATINGS_AGGRESSION = fit_OU_SOCIAL_MATINGS_AGGRESSION,
  OU_SOCIAL_MATINGS_AGGRESSION_INTEGRATION = fit_OU_SOCIAL_MATINGS_AGGRESSION_INTEGRATION
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

###########################
# Platyrrhini: SOCIAL_ORG, MATING SYSTEM, PROP M/F, AGGRESSION
#####################################################
# Platyrrhini
align_filtrado_p = align_filtrado[39:62,]
evolvas_filtrado_p = evolvas_filtrado[39:62,]
mds_p = mds[39:62,]
tree_p = drop.tip(tree, setdiff(tree$tip.label, align_filtrado_p$matings.especies))
dimorfismo_mat_p = dimorfismo_mat[39:62,]

#
df = data.frame(mds_p$especies, mds_p$dados.PARVORDER, mds_p$dados.SOCIAL_ORGANIZATION,
                mds_p$dados.MATING_SYSTEM, mds_p$dados.PROP_MALES_FEMALES, mds_p$dados.DOMINANCE,
                mds_p$dados.AGGRESSION, mds_p$dados.ALL_MALE_GROUPS, mds_p$dados.FURTIVE_COPULATION,
                mds_p$dados.INFANTICIDE, mds_p$dados.MULTILEVEL_SOCIETY, mds_p$dados.N_PAPERS,
                align_filtrado_p$dimor, align_filtrado_p$normas, evolvas_filtrado_p$Size, 
                evolvas_filtrado_p$Integration)

df <- df[match(tree_p$tip.label, df$mds_p.especies), ]

# Modelos
# BM
fit_BM <- BrownianMotionModel(
  phyltree = tree_p,
  mData = dimorfismo_mat_p
)
save(fit_BM, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")

# OU Simple
fit_OU <- ouchModel(
  phyltree = tree_p,
  mData = dimorfismo_mat_p,
  Atype = "Diagonal",        
  estimate.root.state = FALSE
)
save(fit_OU, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_SIMPLE_OU.RData")

# OU SOCIAL ORG
mData_OU_Social <- cbind(dimorfismo_mat_p, df$mds_p.dados.SOCIAL_ORGANIZATION)
colnames(mData_OU_Social)[40] <- "SOCIAL"
rownames(mData_OU_Social) <- df$mds_p.especies
mData_OU_Social <- mData_OU_Social[tree_p$tip.label, ]

fit_OU_SOCIAL <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Social,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_SOCIAL.RData")

# OU MATING SYSTEM
mData_OU_Mating <- cbind(dimorfismo_mat_p, df$mds_p.dados.MATING_SYSTEM)
colnames(mData_OU_Mating)[40] <- "MATING"
rownames(mData_OU_Mating) <- df$mds_p.especies
mData_OU_Mating <- mData_OU_Mating[tree_p$tip.label, ]

fit_OU_MATINGS <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Mating,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_MATINGS, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_MATINGS.RData")

# OU PROP M/F
mData_OU_Prop <- cbind(dimorfismo_mat_p, df$mds_p.dados.PROP_MALES_FEMALES)
colnames(mData_OU_Prop)[40] <- "PROP"
rownames(mData_OU_Prop) <- df$mds_p.especies
mData_OU_Prop <- mData_OU_Prop[tree_p$tip.label, ]

fit_OU_PROP <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Prop,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_PROP, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_PROP.RData")

# OU AGGRESSION
mData_OU_Aggression <- cbind(dimorfismo_mat_p, df$mds_p.dados.AGGRESSION)
colnames(mData_OU_Aggression)[40] <- "AGGRESSION"
rownames(mData_OU_Aggression) <- df$mds_p.especies
mData_OU_Aggression <- mData_OU_Aggression[tree_p$tip.label, ]

fit_OU_AGGRESSION <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Aggression,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_AGGRESSION, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_AGGRESSION.RData")

# OU INTEGRATION
mData_OU_Integration <- cbind(dimorfismo_mat_p, df$evolvas_filtrado_p.Integration)
colnames(mData_OU_Integration)[40] <- "INTEGRATION"
rownames(mData_OU_Integration) <- df$mds_p.especies
mData_OU_Integration <- mData_OU_Integration[tree_p$tip.label, ]

fit_OU_INTEGRATION <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Integration,
  kY = 39,
  predictors = 40, 
  Atype = "Diagonal"
)
save(fit_OU_INTEGRATION, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_INTEGRATION.RData")

# OU SOCIAL + MATINGS
mData_OU_Social_Matings <- cbind(dimorfismo_mat_p, 
                                 df$mds_p.dados.SOCIAL_ORGANIZATION,
                                 df$mds_p.dados.MATING_SYSTEM)
colnames(mData_OU_Social_Matings)[40:41] <- c("SOCIAL", "MATINGS")
rownames(mData_OU_Social_Matings) <- df$mds_p.especies
mData_OU_Social_Matings <- mData_OU_Social_Matings[tree_p$tip.label, ]

fit_OU_SOCIAL_MATINGS <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Social_Matings,
  kY = 39,
  predictors = c(40, 41),
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL_MATINGS, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_SOCIAL_MATINGS.RData")

# OU SOCIAL + MATINGS + PROP
mData_OU_Social_Matings_Prop <- cbind(dimorfismo_mat_p, 
                                 df$mds_p.dados.SOCIAL_ORGANIZATION,
                                 df$mds_p.dados.MATING_SYSTEM,
                                 df$mds_p.dados.PROP_MALES_FEMALES)
colnames(mData_OU_Social_Matings_Prop)[40:42] <- c("SOCIAL", "MATINGS", "PROP")
rownames(mData_OU_Social_Matings_Prop) <- df$mds_p.especies
mData_OU_Social_Matings_Prop <- mData_OU_Social_Matings_Prop[tree_p$tip.label, ]

fit_OU_SOCIAL_MATINGS_PROP <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Social_Matings_Prop,
  kY = 39,
  predictors = c(40, 42),
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL_MATINGS_PROP, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_SOCIAL_MATINGS_PROP.RData")

# OU SOCIAL + MATINGS + PROP + AGGRESSION
mData_OU_Social_Matings_Prop_Aggression <- cbind(dimorfismo_mat_p, 
                                      df$mds_p.dados.SOCIAL_ORGANIZATION,
                                      df$mds_p.dados.MATING_SYSTEM,
                                      df$mds_p.dados.PROP_MALES_FEMALES,
                                      df$mds_p.dados.AGGRESSION)
colnames(mData_OU_Social_Matings_Prop_Aggression)[40:43] <- c("SOCIAL", "MATINGS", "PROP", "AGGRESSION")
rownames(mData_OU_Social_Matings_Prop_Aggression) <- df$mds_p.especies
mData_OU_Social_Matings_Prop_Aggression <- mData_OU_Social_Matings_Prop_Aggression[tree_p$tip.label, ]

fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Social_Matings_Prop_Aggression,
  kY = 39,
  predictors = c(40, 43),
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION.RData")

# OU SOCIAL + MATINGS + PROP + AGGRESSION + INTEGRATION
mData_OU_Social_Matings_Prop_Aggression_Integration <- cbind(dimorfismo_mat_p, 
                                                 df$mds_p.dados.SOCIAL_ORGANIZATION,
                                                 df$mds_p.dados.MATING_SYSTEM,
                                                 df$mds_p.dados.PROP_MALES_FEMALES,
                                                 df$mds_p.dados.AGGRESSION,
                                                 df$evolvas_filtrado_p.Integration)
colnames(mData_OU_Social_Matings_Prop_Aggression_Integration)[40:44] <- c("SOCIAL", "MATINGS", "PROP", "AGGRESSION", "INTEGRATION")
rownames(mData_OU_Social_Matings_Prop_Aggression_Integration) <- df$mds_p.especies
mData_OU_Social_Matings_Prop_Aggression_Integration <- mData_OU_Social_Matings_Prop_Aggression_Integration[tree_p$tip.label, ]

fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_Social_Matings_Prop_Aggression_Integration,
  kY = 39,
  predictors = c(40, 44),
  Atype = "Diagonal"
)
save(fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION.RData")

fits = list(
  BM = fit_BM,
  OU_SIMPLE = fit_OU,
  OU_SOCIAL = fit_OU_SOCIAL,
  OU_MATINGS = fit_OU_MATINGS,
  OU_PROP = fit_OU_PROP,
  OU_AGGRESSION = fit_OU_AGGRESSION,
  OU_INTEGRATION = fit_OU_INTEGRATION,
  OU_SOCIAL_MATINGS = fit_OU_SOCIAL_MATINGS,
  OU_SOCIAL_MATINGS_PROP = fit_OU_SOCIAL_MATINGS_PROP,
  OU_SOCIAL_MATINGS_PROP_AGGRESSION = fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION,
  OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION = fit_OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION
)

save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fits.RData")