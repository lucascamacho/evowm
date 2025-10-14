# setwd
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts")

if(!require(ape)){install.packages("ape"); library(ape)}
if(!require(phytools)){install.packages("phytools"); library(phytools)}
if(!require(geiger)){install.packages("geiger"); library(geiger)}
if(!require(BAMMtools)){install.packages("BAMMtools"); library(BAMMtools)}
if(!require(psych)){install.packages("psych"); library(psych)}
library(readxl)
library(dplyr)

# read and clean data
dados = read_xls("~/Dropbox/Doc/Code/evowm/R/Outputs/NichesPapers_Variables.xls")

cols <- c("Esp$", "DM_geomeanrate", "CVn", "climacv")

dados <- dados %>%
  filter(if_all(all_of(cols), ~ .x != "" & !is.na(.x)))

#Load the phylogeny
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
species = dados$`Esp$`
tree = drop.tip(tree, setdiff(tree$tip.label, species))

to_remove = setdiff(species, tree$tip.label)
dados <- dados %>%
  filter(!`Esp$` %in% to_remove)

#
#pdf("~/Dropbox/Nicho_Dimorfismo_Clima.pdf")
align.vec = setNames(dados$CVn, dados$`Esp$`)
anc <- fastAnc(tree, align.vec)
align.contMap <- contMap(tree, align.vec, fsize = 0.9, lwd = 1, outline = FALSE)
title("CVn", line = -1)  
nodelabels(round(anc, 3), frame = "none", cex = 0.6, adj = c(1.2, -0.4))
#dev.off()

#pdf("~/Dropbox/Plot_2.pdf")
align.vec = setNames(dados$climacv, dados$`Esp$`)
anc <- fastAnc(tree, align.vec)
align.contMap <- contMap(tree, align.vec, fsize = 0.9, lwd = 1, outline = FALSE)
title("climacv", line = -1)  
nodelabels(round(anc, 3), frame = "none", cex = 0.6, adj = c(1.2, -0.4))
#dev.off()

#pdf("~/Dropbox/Plot_3.pdf")
align.vec = setNames(dados$DM_geomeanrate, dados$`Esp$`)
anc <- fastAnc(tree, align.vec)
align.contMap <- contMap(tree, align.vec, fsize = 0.9, lwd = 1, outline = FALSE)
title("DM_geomeanrate", line = -1)  
nodelabels(round(anc, 3), frame = "none", cex = 0.6, adj = c(1.2, -0.4))
#dev.off()

# PICS
dados <- as.data.frame(dados)
rownames(dados) = dados$`Esp$`

#pairs.panels(dados, method = "spearman", hist.col = "#cceeff")

# PIC
dados_ordenado <- dados[tree$tip.label, ]

dados_pic <- as.data.frame(lapply(dados_ordenado[,c(13,14,15,16,17,18,19,20,21,22,23,24,30)], function(x) pic(x, tree)))

#pairs.panels(dados_pic, method = "pearson", hist.col = "#cceeff")
corr.test(dados_pic, method = "pearson", adjust = "none")

x = cor(dados_pic)
write.xlsx(x, "~/Dropbox/Doc/Code/evowm/R/Outputs/PIC_Correlations_Monkeys_ArgentinaCongress.xlsx")

##### BATS
dados = read_xls("~/Dropbox/Doc/Gabriel_Bats/bats.xls")

#
filename = "~/Dropbox/Doc/Gabriel_Bats/filostTree_nova.nex"
tree = ape::read.nexus(filename)
species = dados$`Esp$`
tree = drop.tip(tree, setdiff(tree$tip.label, species))

bats_data = dados[,c(1,9,10,11,16,17,20,21,25)]

pairs.panels(bats_data[,2:9], method = "pearson", hist.col = "#cceeff")

# PIC
bats_data <- as.data.frame(bats_data)
rownames(bats_data) = dados$`Esp$`
dados_ordenado <- bats_data[tree$tip.label, ]

#
dados_pic <- as.data.frame(lapply(dados_ordenado[,2:9], function(x) pic(x, tree)))

#pairs.panels(dados_pic, method = "pearson", hist.col = "#cceeff")
corr.test(dados_pic, method = "pearson", adjust = "none")

x = cor(dados_pic)
write.xlsx(x, "~/Dropbox/Doc/Code/evowm/R/Outputs/PIC_Correlations_Bats_ArgentinaCongress.xlsx")
