# load packages and functions
setwd("~/Dropbox/Doc/Data/")

if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(ape)){install.packages("ape"); library(ape)}
if(!require(phytools)){install.packages("phytools"); library(phytools)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# load mating systems data 
matings = read.table("~/Dropbox/Doc/Data/wos_mating_systems/primates_mating_systems.csv",
                     sep = ",", header = TRUE)
# get "interested" species
species_tree = vector()
for(i in 1:nrow(matings)){
  if(is.na(matings$SUBSPECIES[i])){
    species_tree[i] = paste(matings$GENUS[i], matings$SPECIES[i], sep="_")}
  else{
    species_tree[i] = paste(matings$GENUS[i], matings$SPECIES[i], matings$SUBSPECIES[i], sep="_")
  }
}

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
cat.tree = keep.tip(tree, species_tree)

# plot tree
plot(cat.tree)
nodelabels(frame = "n", cex = 0.8, col = "blue")

# open dimorphism estimations
machos = read.table("~/Dropbox/Doc/Output/Avgs_Male.txt", header = TRUE)
femeas = read.table("~/Dropbox/Doc/Output/Avgs_Female.txt", header = TRUE)
dimorp = log(machos[,2:40]) - log(femeas[,2:40])
dimorp = cbind(machos$species_tree, dimorp)
colnames(dimorp)[1] = "species_tree"

# estimate ancestral dimorphisms in ML
ancs = sapply(dimorp[,2:40], function(x) ace(x, cat.tree)$ace)

ancs = cbind(seq(72,141,1), ancs)
colnames(ancs)[1] = "species_tree"
ancs_2 = rbind(dimorp, ancs)

# transform all dimorphisms (actual and ancestral) in a table
lista_medias = dlply(ancs_2, .(species_tree), numcolwise(identity))

# read V/CV matrices
setwd("~/Dropbox/Doc/Output/log_vcv")

# read all vcv matrices
temp = list.files(pattern="*.txt")
vcv = lapply(temp, read.table, row.names = 1, header = TRUE)
names(vcv)  = gsub(".txt", replacement= "", temp)

# estimate V/CV matrices from ancestral states
all_cov_matrices = PhyloW(cat.tree, vcv)
attributes(all_cov_matrices)$split_labels = attributes(all_cov_matrices)$names

# calculate PCs from ancestral and actual V/CV matrices
all_pc1 = llply(all_cov_matrices, function(x) eigen(x)$vectors[,1])
all_pc2 = llply(all_cov_matrices, function(x) eigen(x)$vectors[,2])
all_pc3 = llply(all_cov_matrices, function(x) eigen(x)$vectors[,3])
all_pc4 = llply(all_cov_matrices, function(x) eigen(x)$vectors[,4])
all_pc5 = llply(all_cov_matrices, function(x) eigen(x)$vectors[,5])
all_pc6 = llply(all_cov_matrices, function(x) eigen(x)$vectors[,6])

checar = data.frame()
# check nodes PCs, Dimorphisms, Isometric 
for(i in 72:141){
  position_d = which(names(lista_medias) == i)
  
  dimorps = lista_medias[position_d]
  nor = norma(unlist(dimorps))
  
  position_1 = which(names(all_pc1) == i)
  position_2 = which(names(all_pc2) == i)
  position_3 = which(names(all_pc3) == i)
  position_4 = which(names(all_pc4) == i)
  
  dimorp_PC1 = abs(corVector(unlist(all_pc1[position_1]), unlist(dimorps)))
  dimorp_PC2 = abs(corVector(unlist(all_pc2[position_2]), unlist(dimorps)))
  dimorp_PC3 = abs(corVector(unlist(all_pc3[position_3]), unlist(dimorps)))
  dimorp_PC4 = abs(corVector(unlist(all_pc4[position_4]), unlist(dimorps)))
  
  iso_PC1 = abs(corVector(unlist(all_pc1[position_1]), rep(0.160128154, 39)))
  iso_PC2 = abs(corVector(unlist(all_pc2[position_2]), rep(0.160128154, 39)))
  iso_PC3 = abs(corVector(unlist(all_pc3[position_3]), rep(0.160128154, 39)))
  iso_PC4 = abs(corVector(unlist(all_pc4[position_4]), rep(0.160128154, 39)))
  
  Dimorp_Norms = nor
  
  lixo = c(i, dimorp_PC1, dimorp_PC2, dimorp_PC3, dimorp_PC4, iso_PC1, iso_PC2, iso_PC3, iso_PC4, Dimorp_Norms)
  checar = rbind(checar, lixo)
  
  # temporario
  
  
}

names(checar) = c("node", "dimorp_PC1", "dimorp_PC2", "dimorp_PC3", "dimorp_PC4", "iso_PC1", "iso_PC2", "iso_PC3", "iso_PC4", "Dimorp_Norms")

write.table(checar, file = "~/Dropbox/Doc/Output/Log_Ancestrals_Pcs_Dimorp_Norma.txt", row.names = FALSE, dec = ".", sep = '\t', quote = FALSE)

#
# calculate Delta_Z
node_names = c(cat.tree$tip.label, 72:141)
node_numbers = c(1:71, 73:141)
delta_Zs = vector("list", length(node_numbers))

for(i in 1:length(node_numbers)){
  node = node_numbers[i]
  current_node = node_names[node]
  node_mean = lista_medias[[current_node]]
  ancestral_node = node_names[cat.tree$edge[cat.tree$edge[,2] == node, 1]] #
  ancestral_mean = lista_medias[[ancestral_node]]
  delta_Zs[[i]] = node_mean - ancestral_mean
}

names(delta_Zs) = node_names[-72]
lista_deltaZs = ldply(delta_Zs)

# save table
#write.table(lista_deltaZs, file = "~/Dropbox/Doc/Output/Delta_Zs.txt", row.names = FALSE, 
#            dec = ".", sep = '\t', quote = FALSE)

# cor pcs e delta Zs
cor_deltz_pcs = data.frame()
for(i in 1:nrow(lista_deltaZs)){
  # choose species Delta_Z with the ancestral
  name_node = lista_deltaZs$.id[i]
  
  # Cor Vector between Delta_Z and PCs
  pcum = abs(corVector(lista_deltaZs[i,2:40], unlist(all_pc1[grepl(name_node, names(all_pc1))])))
  pdois = abs(corVector(lista_deltaZs[i,2:40], unlist(all_pc2[grepl(name_node, names(all_pc2))])))
  ptres = abs(corVector(lista_deltaZs[i,2:40], unlist(all_pc3[grepl(name_node, names(all_pc3))])))
  pquatro = abs(corVector(lista_deltaZs[i,2:40], unlist(all_pc4[grepl(name_node, names(all_pc4))])))
  
  # Cor Vector between PCs and Isometric vector
  cor_iso_1 = abs(corVector(unlist(all_pc1[grepl(name_node, names(all_pc1))]), rep(0.160128154, 39)))
  cor_iso_2 = abs(corVector(unlist(all_pc2[grepl(name_node, names(all_pc2))]), rep(0.160128154, 39)))
  
  # get all data together
  data = c(pcum, pdois, ptres, pquatro, cor_iso_1, cor_iso_2)
  cor_deltz_pcs = rbind(cor_deltz_pcs, data)
}

# organize and rename collum names
cor_deltz_pcs = data.frame(names(delta_Zs), cor_deltz_pcs)
colnames(cor_deltz_pcs) = c("Nodes", "Cor_PC1_DeltaZ", "Cor_PC2_DeltaZ", "Cor_PC3_DeltaZ", "Cor_PC4_DeltaZ",
                            "Cor_Iso_PC1", "Cor_Iso_PC2")

# save correlation of PCs and Delta_Z table
#write.table(cor_deltz_pcs, file = "~/Dropbox/Doc/Output/Corr_DeltaZ_PCs.txt", row.names = FALSE, 
#            dec = ".", sep = '\t', quote = FALSE)
