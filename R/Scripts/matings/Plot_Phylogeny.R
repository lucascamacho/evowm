setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

#
library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
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

#################################################################################
## TREE
#################################################################################
tree = tree_genus 
dimorfismo = data.frame(genero = align_f$matings.genus, 
                        Parvorder = align_f$matings.dados.PARVORDER, 
                        normas = align_f$normas)

mycols = c("#1b9e77", "#7570b3")

# vetores com os gêneros de cada grupo
cats <- dimorfismo$genero[dimorfismo$Parvorder == "Catarrhini"]
plats <- dimorfismo$genero[dimorfismo$Parvorder == "Platyrrhini"]

# nós MRCA
node_cat <- getMRCA(tree, cats)
node_plat <- getMRCA(tree, plats)

# conferir nomes
setdiff(tree$tip.label, dimorfismo$genero)

# ordenar o df para bater com a árvore (opcional mas ajuda)
dimorfismo <- dimorfismo %>%
  filter(genero %in% tree$tip.label)

p <- ggtree(tree, layout = "circular", linewidth = 1.2, color = "black")

p$data <- p$data %>%
  left_join(dimorfismo, by = c("label" = "genero"))

mycols <- c("Catarrhini" = "#1b9e77",
            "Platyrrhini" = "#7570b3")

p_final <- p +
  # colorir tips por parvorder
  geom_tippoint(aes(color = Parvorder), size = 2) +
  
  # nomes dos gêneros
  geom_tiplab(size = 4, offset = 2.5, align = TRUE) +
  
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

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_Phylogeny.png",
       plot = p_final,
       width = 10,
       height = 10,
       units = "in",
       dpi = 600)
