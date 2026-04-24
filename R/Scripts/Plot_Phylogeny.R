setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

#
library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
library(dplyr)

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

#################################################################################
## TREE
#################################################################################
dimorfismo = data.frame(genero = mds_f$dados.GENUS,
                        species = align_f$matings.especies, 
                        Parvorder = mds_f$dados.PARVORDER, 
                        normas = align_f$normas)

mycols = c("#1b9e77", "#7570b3")

# vetores com os gêneros de cada grupo
cats <- dimorfismo$species[dimorfismo$Parvorder == "Catarrhini"]
plats <- dimorfismo$species[dimorfismo$Parvorder == "Platyrrhini"]

# nós MRCA
node_cat <- getMRCA(tree, cats)
node_plat <- getMRCA(tree, plats)

# conferir nomes
setdiff(tree$tip.label, dimorfismo$species)

# ordenar o df para bater com a árvore (opcional mas ajuda)
dimorfismo <- dimorfismo %>%
  filter(species %in% tree$tip.label)

p <- ggtree(tree, layout = "circular", linewidth = 1.2, color = "black")

p$data <- p$data %>%
  left_join(dimorfismo, by = c("label" = "species"))

mycols <- c("Catarrhini" = "#1b9e77",
            "Platyrrhini" = "#7570b3")

p_final <- p +
  # colorir tips por parvorder
  geom_tippoint(aes(color = Parvorder), size = 2) +
  
  # nomes dos gêneros
  geom_tiplab(size = 3, offset = 2.5, align = TRUE) +
  
  geom_point(aes(x = x + 1, size = normas),
             color = "red",
             na.rm = TRUE) +
  
  scale_color_manual(values = mycols) +
  
  scale_size_continuous(name = "Sexual Dimorphism",
                        range = c(2, 6)) +
  # highlight dos clados
  geom_hilight(node = node_cat, fill = "#1b9e77", alpha = 0.4) +
  geom_hilight(node = node_plat, fill = "#7570b3", alpha = 0.4) +
  
  theme(legend.position = "right",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(0.8, "cm"))

p_final

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_Phylogeny.png",
       plot = p_final,
       width = 14,
       height = 14,
       units = "in",
       dpi = 600)
