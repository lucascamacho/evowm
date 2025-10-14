# Average and PCs of ancestral species
# set WD
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts")

# load functions
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(ape)){install.packages("ape"); library(ape)}
if(!require(phytools)){install.packages("phytools"); library(phytools)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}

# read all VCV matrices
setwd("~/Dropbox/Doc/Data/vcv/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement = "", temp)
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts")

# get species
matings = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")

# Define species based on match between vcv and matings
commom <- intersect(names(vcv), matings$especies)

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
species = matings$especies
tree = drop.tip(tree, setdiff(tree$tip.label, commom))

# get species averages for males and females
averages = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

# calculate male and female averages from phylogeny
avg_males = data.frame()
avg_females = data.frame()
for(i in 1:length(tree$tip.label)){
  name_species = tree$tip.label[i]
  index = which(names(averages$ByTrait_Averages) == name_species)
  med_machos = averages$ByTrait_Averages[[index]]$Machos
  med_femeas = averages$ByTrait_Averages[[index]]$Fêmeas
  
  avg_males = rbind(avg_males, med_machos)
  avg_females = rbind(avg_females, med_femeas)
}

colnames(avg_males) = names(averages$ByTrait_Averages$Cercocebus_agilis$Machos)
avg_males = cbind(tree$tip.label, avg_males)
colnames(avg_males)[1] = "Species"

colnames(avg_females) = names(averages$ByTrait_Averages$Cercocebus_agilis$Fêmeas)
avg_females = cbind(tree$tip.label, avg_females)
colnames(avg_females)[1] = "Species"

# ACE MALES
rownames(avg_males) <- avg_males[,1]
# 
anc_states_male <- sapply(avg_males[, 2:40], function(x) {
  names(x) <- rownames(avg_males)
  fastAnc(tree, x, vars = FALSE, CI = FALSE)  # já é um vetor
})

# ACE FEMALES
rownames(avg_females) <- avg_females[,1]
# 
anc_states_female <- sapply(avg_females[, 2:40], function(x) {
  names(x) <- rownames(avg_females)
  fastAnc(tree, x, vars = FALSE, CI = FALSE)  # já é um vetor
})

# single list
bytrait_averages = list(anc_states_male, anc_states_female)

# calculate dimorphism of species from phylogeny
dimor = vector()
for(i in 1:length(tree$tip.label)){
  name_species = tree$tip.label[i]
  sp = averages$Averages[[name_species]]
  dimor[i] = unlist(sp[1]) - unlist(sp[2])
}

# ACE DIMOR
x = setNames(dimor, avg_females[,1]) 
dimorps_bytrait = fastAnc(tree, x, vars = FALSE, CI = FALSE)  # já é um vetor

# estimate ancestral dimorphism by trait in ML
dimorps_bytrait = data.frame()
for(i in 1:length(tree$tip.label)){
  name_species = tree$tip.label[i]
  sp = averages$ByTrait_Averages[[name_species]]
  dimor_bytrait = unlist(sp[1]) - unlist(sp[2])
  dimorps_bytrait = rbind(dimorps_bytrait, dimor_bytrait)
}

# rename columns and insert species names
colnames(dimorps_bytrait) = names(averages$ByTrait_Averages$Cercocebus_agilis$Machos)
dimorps_bytrait = cbind(tree$tip.label, dimorps_bytrait)
colnames(dimorps_bytrait)[1] = "Species"

# ACE DIMOR BY TRAIT
rownames(dimorps_bytrait) <- avg_females[,1]
# 
bytrait_dimor_ancs <- sapply(dimorps_bytrait[, 2:40], function(x) {
  names(x) <- rownames(avg_females)
  fastAnc(tree, x, vars = FALSE, CI = FALSE)  # já é um vetor
})

# estimate VCV and PCs for ancestral
index = setdiff(names(vcv), species)
vcv2 = vcv[!names(vcv) %in% index]
all_cov_matrices = PhyloW(tree, vcv2)

# PCs of the ancestral species
pcs = list()
vals = list()
diags = list()
for(i in 63:123){
  sp = which(names(all_cov_matrices) == i)
  covar = all_cov_matrices[[sp]]
  
  # get 8 PCs and eigenvalues
  pc_1 = eigen(covar)$vectors[,1]
  vals_1 = eigen(covar)$values[1]
  pc_2 = eigen(covar)$vectors[,2]
  vals_2 = eigen(covar)$values[2]
  pc_3 = eigen(covar)$vectors[,3]
  vals_3 = eigen(covar)$values[3]
  pc_4 = eigen(covar)$vectors[,4]
  vals_4 = eigen(covar)$values[4]
  pc_5 = eigen(covar)$vectors[,5]
  vals_5 = eigen(covar)$values[5]
  pc_6 = eigen(covar)$vectors[,6]
  vals_6 = eigen(covar)$values[6]
  pc_7 = eigen(covar)$vectors[,7]
  vals_7 = eigen(covar)$values[7]
  pc_8 = eigen(covar)$vectors[,8]
  vals_8 = eigen(covar)$values[8]
  
  # diag
  v = diag(as.matrix(covar))
  
  pcs[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6, pc_7, pc_8)
  vals[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6, vals_7, vals_8)
  diags[[i]] = v
}

# make list and rename it
ancestral_dimorp_pcs = list(seq(63, 123, 1), 
                            bytrait_averages, 
                            rowMeans(bytrait_averages[[1]] - bytrait_averages[[2]]), 
                            bytrait_dimor_ancs, 
                            rev(all_cov_matrices[63:123]),  
                            pcs, 
                            vals, 
                            diags)
names(ancestral_dimorp_pcs) = c("Ancestrals", "ByTrait_Averages", "Dimorphism","ByTrait_Dimorphism", "VCV", "PCs", "Eigenvalues", "Diagonals")

# name lists and remove empty lists
names(ancestral_dimorp_pcs$ByTrait_Averages) = c("Machos", "Fêmeas")

names(ancestral_dimorp_pcs$PCs) = seq(1, 123, 1)
names(ancestral_dimorp_pcs$Eigenvalues) = seq(1, 123, 1)
names(ancestral_dimorp_pcs$Diagonals) = seq(1, 123, 1)

ancestral_dimorp_pcs$PCs = ancestral_dimorp_pcs$PCs[-(1:62)]
ancestral_dimorp_pcs$Eigenvalues = ancestral_dimorp_pcs$Eigenvalues[-(1:62)]
ancestral_dimorp_pcs$Diagonals = ancestral_dimorp_pcs$Diagonals[-(1:62)]

for (i in 1:61) {
  names(ancestral_dimorp_pcs$PCs[[i]]) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  names(ancestral_dimorp_pcs$Eigenvalues[[i]]) <- c("Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6", "Lambda7", "Lambda8")
}

# save final list of results
saveRDS(ancestral_dimorp_pcs, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/ancestrals_averages_PCS_autovalues_primates.RDS")
