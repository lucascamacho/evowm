setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(mvSLOUCH)

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

align_filtrado <- align[align$matings.especies %in% mds$especies, ]
evolvas_filtrado <- evolvas[evolvas$species %in% mds$especies, ]

species62 = align_filtrado$medidas.Species

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
align_filtrado_c = align_filtrado[1:38,]
evolvas_filtrado_c = evolvas_filtrado[1:38,]
mds_c = mds[1:38,]
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_filtrado_c$matings.especies))

#
df = data.frame(mds_c$especies, mds_c$scores.fit....1., mds_c$scores.fit....2., mds_c$dados.PARVORDER,
                align_filtrado_c$dimor, align_filtrado_c$normas, align_filtrado_c$align_1, align_filtrado_c$align_2,
                align_filtrado_c$align_3, align_filtrado_c$align_4, align_filtrado_c$align_5, align_filtrado_c$align_6,
                evolvas_filtrado_c$Dimorfism, evolvas_filtrado_c$Integration_A, evolvas_filtrado_c$Size, 
                evolvas_filtrado_c$Evolvability, evolvas_filtrado_c$Conditional_Evolvability, 
                evolvas_filtrado_c$Average_Evolvability)

df <- df[match(tree_c$tip.label, df$mds_c.especies), ]

# Modelos
# BM
mData_BM = matrix(df$align_filtrado_c.dimor, ncol = 1)
colnames(mData_BM) <- "SD"
rownames(mData_BM) <- df$mds_c.especies

fit_BM <- BrownianMotionModel(
  phyltree = tree_c,
  mData = mData_BM
)
save(fit_BM, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_BM.RData")

# OU Simple
mData_OU <- matrix(df$align_filtrado_c.dimor, ncol = 1)
rownames(mData_OU) <- df$mds_c.especies
mData_OU <- mData_OU[tree_c$tip.label, , drop = FALSE]
colnames(mData_OU) <- "SD"

fit_OU <- ouchModel(
  phyltree = tree_c,
  mData = mData_OU,
  Atype = "Diagonal",        
  estimate.root.state = FALSE
)

save(fit_OU, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_ouch_align1.RData")

# OU NMDS1 - 
mData_OU_MDS1 <- cbind(df$align_filtrado_c.dimor, df$mds_c.scores.fit....1.)
colnames(mData_OU_MDS1) <- c("SD", "NMDS1")
rownames(mData_OU_MDS1) <- df$mds_c.especies
mData_OU_MDS1 <- mData_OU_MDS1[tree_c$tip.label, ]

fit_OU_MDS1 <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_MDS1,
  kY = 1,
  predictors = 2,  # coluna NMDS1 como preditor
  Atype = "Diagonal"
)
save(fit_OU_MDS1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_MDS1.RData")

# OU NMDS2
mData_OU_MDS2 <- cbind(df$align_filtrado_c.dimor,
                       df$mds_c.scores.fit....2.)
colnames(mData_OU_MDS2) <- c("SD", "NMDS2")
rownames(mData_OU_MDS2) <- df$mds_c.especies
mData_OU_MDS2 <- mData_OU_MDS2[tree_c$tip.label, ]
fit_OU_MDS2 <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_MDS2,
  kY = 1,
  predictors = 2,
  Atype = "Diagonal"
)
save(fit_OU_MDS2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_MDS2.RData")

# OU NMDS1 + NMDS2
mData_OU_MDS12 <-  cbind(df$align_filtrado_c.dimor, 
                         df$mds_c.scores.fit....1.,
                         df$mds_c.scores.fit....2.)
colnames(mData_OU_MDS12) <- c("SD","NMDS1", "NMDS2")
rownames(mData_OU_MDS12) <- df$mds_c.especies
mData_OU_MDS12 <- mData_OU_MDS12[tree_c$tip.label, ]
fit_OU_MDS12 <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_MDS12,
  kY = 1,
  predictors = 2:3,
  Atype = "Diagonal")
save(fit_OU_MDS12, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_MDS12.RData")

# OU NMDS1 + NMDS2 + Size
mData_OU_MDS12_size <- cbind(df$align_filtrado_c.dimor,
                             df$mds_c.scores.fit....1.,
                             df$mds_c.scores.fit....2.,
                             df$evolvas_filtrado_c.Size)
colnames(mData_OU_MDS12_size) <- c("SD", "NMDS1", "NMDS2", "Size")
rownames(mData_OU_MDS12_size) <- df$mds_c.especies
mData_OU_MDS12_size <- mData_OU_MDS12_size[tree_c$tip.label, ]
fit_OU_MDS12_size <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_MDS12_size,
  kY = 1,
  predictors = 2:4,
  Atype = "Diagonal"
)
save(fit_OU_MDS12_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_MDS12_size.RData")

# OU NMDS1 + NMDS2 + Size + Cond Evolvability
mData_OU_MDS12_size_cond_evolv <- cbind(df$align_filtrado_c.dimor, 
                             df$mds_c.scores.fit....1.,
                             df$mds_c.scores.fit....2.,
                             df$evolvas_filtrado_c.Size,
                             df$evolvas_filtrado_c.Conditional_Evolvability)
colnames(mData_OU_MDS12_size_cond_evolv) <- c("SD","NMDS1", "NMDS2", "Size", "Cond_Evolv")
rownames(mData_OU_MDS12_size_cond_evolv) <- df$mds_c.especies
mData_OU_MDS12_size_cond_evolv <- mData_OU_MDS12_size_cond_evolv[tree_c$tip.label, ]
fit_OU_OU_MDS12_size_cond_evolv <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_MDS12_size_cond_evolv,
  kY = 1,
  predictors = 2:5,
  Atype = "Diagonal"
)
save(fit_OU_MDS12_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_OU_MDS12_size_cond_evolv.RData")

# OU Full
mData_OU_full <- cbind(df$align_filtrado_c.dimor,
                       df$mds_c.scores.fit....1.,
                       df$mds_c.scores.fit....2.,
                       df$evolvas_filtrado_c.Size,
                       df$evolvas_filtrado_c.Conditional_Evolvability,
                       df$evolvas_filtrado_c.Integration_A)
colnames(mData_OU_full) <- c("SD","NMDS1", "NMDS2", "Size", "Cond_Evolv", "Integration")
rownames(mData_OU_full) <- df$mds_c.especies
mData_OU_full <- mData_OU_full[tree_c$tip.label, ]
fit_OU_full <- mvslouchModel(
  phyltree = tree_c,
  mData = mData_OU_full,
  kY = 1,
  predictors = 2:6,
  Atype = "Diagonal"
)
save(fit_OU_full, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fit_OU_full.RData")

fits <- list(
  BM = fit_BM,
  Simple_OU = fit_OU,
  OU_MDS1 = fit_OU_MDS1,
  OU_MDS2 = fit_OU_MDS2,
  OU_MDS12 = fit_OU_MDS12,
  OU_MDS12_size = fit_OU_MDS12_size,
  OU_MDS12_size_cond_evolvability = fit_OU_OU_MDS12_size_cond_evolv,
  OU_full = fit_OU_full
)
save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_mvslouch_fits.RData")

# Platyrrhini
align_filtrado_p = align_filtrado[39:62,]
evolvas_filtrado_p = evolvas_filtrado[39:62,]
mds_p = mds[39:62,]
tree_p = drop.tip(tree, setdiff(tree$tip.label, align_filtrado_p$matings.especies))

#
df = data.frame(mds_p$especies, mds_p$scores.fit....1., mds_p$scores.fit....2., mds_p$dados.PARVORDER,
                align_filtrado_p$dimor, align_filtrado_p$normas, align_filtrado_p$align_1, align_filtrado_p$align_2,
                align_filtrado_p$align_3, align_filtrado_p$align_4, align_filtrado_p$align_5, align_filtrado_p$align_6,
                evolvas_filtrado_p$Dimorfism, evolvas_filtrado_p$Integration_A, evolvas_filtrado_p$Size, 
                evolvas_filtrado_p$Evolvability, evolvas_filtrado_p$Conditional_Evolvability, 
                evolvas_filtrado_p$Average_Evolvability)

df <- df[match(tree_p$tip.label, df$mds_p.especies), ]

# Modelos
# BM
mData_BM = matrix(df$align_filtrado_p.dimor, ncol = 1)
colnames(mData_BM) <- "SD"
rownames(mData_BM) <- df$mds_p.especies

fit_BM <- BrownianMotionModel(
  phyltree = tree_p,
  mData = mData_BM
)

save(fit_BM, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_BM.RData")

# OU Simple
mData_OU <- matrix(df$align_filtrado_p.dimor, ncol = 1)
rownames(mData_OU) <- df$mds_p.especies
mData_OU <- mData_OU[tree_p$tip.label, , drop = FALSE]
colnames(mData_OU) <- "SD"

fit_OU <- ouchModel(
  phyltree = tree_p,
  mData = mData_OU,
  Atype = "Diagonal",        
  estimate.root.state = FALSE
)

save(fit_OU, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_ouch_align1.RData")

# OU NMDS1 - 
mData_OU_MDS1 <- cbind(df$align_filtrado_p.dimor, df$mds_p.scores.fit....1.)
colnames(mData_OU_MDS1) <- c("SD", "NMDS1")
rownames(mData_OU_MDS1) <- df$mds_p.especies
mData_OU_MDS1 <- mData_OU_MDS1[tree_p$tip.label, ]

fit_OU_MDS1 <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_MDS1,
  kY = 1,
  predictors = 2,  # coluna NMDS1 como preditor
  Atype = "Diagonal"
)
save(fit_OU_MDS1, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_MDS1.RData")

# OU NMDS2
mData_OU_MDS2 <- cbind(df$align_filtrado_p.dimor,
                       df$mds_p.scores.fit....2.)
colnames(mData_OU_MDS2) <- c("SD", "NMDS2")
rownames(mData_OU_MDS2) <- df$mds_p.especies
mData_OU_MDS2 <- mData_OU_MDS2[tree_p$tip.label, ]
fit_OU_MDS2 <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_MDS2,
  kY = 1,
  predictors = 2,
  Atype = "Diagonal"
)
save(fit_OU_MDS2, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_MDS2.RData")

# OU NMDS1 + NMDS2
mData_OU_MDS12 <-  cbind(df$align_filtrado_p.dimor, 
                         df$mds_p.scores.fit....1.,
                         df$mds_p.scores.fit....2.)
colnames(mData_OU_MDS12) <- c("SD","NMDS1", "NMDS2")
rownames(mData_OU_MDS12) <- df$mds_p.especies
mData_OU_MDS12 <- mData_OU_MDS12[tree_p$tip.label, ]
fit_OU_MDS12 <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_MDS12,
  kY = 1,
  predictors = 2:3,
  Atype = "Diagonal")
save(fit_OU_MDS12, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_MDS12.RData")

# OU NMDS1 + NMDS2 + Size
mData_OU_MDS12_size <- cbind(df$align_filtrado_p.dimor,
                             df$mds_p.scores.fit....1.,
                             df$mds_p.scores.fit....2.,
                             df$evolvas_filtrado_p.Size)
colnames(mData_OU_MDS12_size) <- c("SD", "NMDS1", "NMDS2", "Size")
rownames(mData_OU_MDS12_size) <- df$mds_p.especies
mData_OU_MDS12_size <- mData_OU_MDS12_size[tree_p$tip.label, ]
fit_OU_MDS12_size <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_MDS12_size,
  kY = 1,
  predictors = 2:4,
  Atype = "Diagonal"
)
save(fit_OU_MDS12_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_MDS12_size.RData")

# OU NMDS1 + NMDS2 + Size + Cond Evolvability
mData_OU_MDS12_size_cond_evolv <- cbind(df$align_filtrado_p.dimor, 
                                        df$mds_p.scores.fit....1.,
                                        df$mds_p.scores.fit....2.,
                                        df$evolvas_filtrado_p.Size,
                                        df$evolvas_filtrado_p.Conditional_Evolvability)
colnames(mData_OU_MDS12_size_cond_evolv) <- c("SD","NMDS1", "NMDS2", "Size", "Cond_Evolv")
rownames(mData_OU_MDS12_size_cond_evolv) <- df$mds_p.especies
mData_OU_MDS12_size_cond_evolv <- mData_OU_MDS12_size_cond_evolv[tree_p$tip.label, ]
fit_OU_OU_MDS12_size_cond_evolv <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_MDS12_size_cond_evolv,
  kY = 1,
  predictors = 2:5,
  Atype = "Diagonal"
)
save(fit_OU_MDS12_size, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_OU_MDS12_size_cond_evolv.RData")

# OU Full
mData_OU_full <- cbind(df$align_filtrado_p.dimor,
                       df$mds_p.scores.fit....1.,
                       df$mds_p.scores.fit....2.,
                       df$evolvas_filtrado_p.Size,
                       df$evolvas_filtrado_p.Conditional_Evolvability,
                       df$evolvas_filtrado_p.Integration_A)
colnames(mData_OU_full) <- c("SD","NMDS1", "NMDS2", "Size", "Cond_Evolv", "Integration")
rownames(mData_OU_full) <- df$mds_p.especies
mData_OU_full <- mData_OU_full[tree_p$tip.label, ]
fit_OU_full <- mvslouchModel(
  phyltree = tree_p,
  mData = mData_OU_full,
  kY = 1,
  predictors = 2:6,
  Atype = "Diagonal"
)
save(fit_OU_full, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fit_OU_full.RData")

fits <- list(
  BM = fit_BM,
  Simple_OU = fit_OU,
  OU_MDS1 = fit_OU_MDS1,
  OU_MDS2 = fit_OU_MDS2,
  OU_MDS12 = fit_OU_MDS12,
  OU_MDS12_size = fit_OU_MDS12_size,
  OU_MDS12_size_cond_evolvability = fit_OU_OU_MDS12_size_cond_evolv,
  OU_full = fit_OU_full
)
save(fits, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_mvslouch_fits.RData")
