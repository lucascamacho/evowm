# SLOUCH POR GENERO
# CATARRHINI, PLATYRRHINI E HAPLORRHYNI

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

library(OUwie)
library(ape)
library(dplyr)

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")
nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
  nicho$genus,
  align$matings.genus,
  evolvas$genus,
  mds$genus
))

nicho_f   <- nicho[nicho$genus %in% sp_comuns, ]
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

media_genero <- nicho %>%
  group_by(genus) %>%
  summarise(
    climacv_mean = mean(climacv, na.rm = TRUE),
    CVn_mean = mean(CVn, na.rm = TRUE)
  )

nicho_f <- media_genero[match(tree_genus$tip.label, media_genero$genus), ]
evolvas_f <- evolvas_f[match(tree_genus$tip.label, evolvas_f$genus), ]
align_f <- align_f[match(tree_genus$tip.label, align_f$matings.genus), ]
mds_f <- mds_f[match(tree_genus$tip.label, mds_f$genus), ]

################################################################################
# A PARTIR DAQUI
# Catarrhini 1:16 com nicho, 1:21 sem nicho
align_c = align_f[1:16,]
evolvas_c = evolvas_f[1:16,]
mds_c = mds_f[1:16,]
tree_c = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_c$matings.genus))
nicho_c = nicho_f[1:16,]

#
df = data.frame(mds_c$genus, 
                mds_c$dados.PARVORDER, 
                mds_c$dados.SOCIAL_ORGANIZATION,
                mds_c$dados.MATING_SYSTEM, 
                mds_c$dados.AGGRESSION, evolvas_c$Dimorfism, 
                evolvas_c$Integration, 
                evolvas_c$Size, 
                mds_c$vegan..scores.fit..sites...1., 
                mds_c$vegan..scores.fit..sites...2., 
                align_c$normas,
                nicho_c$climacv_mean, 
                nicho_c$CVn_mean)

df <- df[match(tree_c$tip.label, df$mds_c.genus), ]

df_ou <- data.frame(
  genus = factor(df$mds_c.genus, levels = tree_c$tip.label),
  integration = as.factor(df$evolvas_c.Integration),
  response = align_c$normas
)

df_ou <- df_ou[match(tree_c$tip.label, df_ou$genus), ]

fit_ou_integration = OUwie(phy = tree_c,
                       data = df_ou,
                       models)


#################################################################################
# Platyrrhini 22:37 sem nicho, 17:28 com nicho
align_p = align_f[22:37,]
evolvas_p = evolvas_f[22:37,]
mds_p = mds_f[22:37,]
tree_p = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_p$matings.genus))
nicho_p = nicho_f[17:28,]





#################################################################################
# Haplorrhini
df = data.frame(mds_f$genus, mds_f$dados.PARVORDER, mds_f$dados.SOCIAL_ORGANIZATION,
                mds_f$dados.MATING_SYSTEM, mds_f$dados.AGGRESSION, evolvas_f$Dimorfism, 
                evolvas_f$Integration, 
                evolvas_f$Size, 
                mds_f$vegan..scores.fit..sites...1., 
                mds_f$vegan..scores.fit..sites...2., align_f$normas
                , nicho_f$climacv_mean, nicho_f$CVn_mean)

df <- df[match(tree_genus$tip.label, df$mds_f.genus), ]
