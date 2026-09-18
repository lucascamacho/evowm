# Average and PCs of ancestral species
# log and no-log averages and PCs
# set WD

# load functions
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(ape)){install.packages("ape"); library(ape)}
if(!require(phytools)){install.packages("phytools"); library(phytools)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}

# read species names
setwd("~/Dropbox/Doc/Output/nolog_vcv")
temp = list.files(pattern = "*.txt")
names(temp)  = gsub(".txt", replacement= "", temp)

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)

# remove species which are not in the phylogeny
to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), names(temp))
temp = temp[-index]
cat.tree = keep.tip(tree, names(temp))

# get species averages for males and females
averages = readRDS("~/Dropbox/Novo_Output/Averages_PCS_Extant_Species.RData")

# calculate dimorphism
dimor_nolog = vector()
dimor_log = vector()

for(i in 1:length(averages$Species)){
  sp = averages$NoLog_Averages[[i]]
  dimor_nolog[i] = unlist(sp[1]) - unlist(sp[2])
  
  sp = averages$Log_Averages[[i]]
  dimor_log[i] = unlist(sp[1]) - unlist(sp[2])
}

# estimate ancestral dimorphism in ML
dimor_ancs_log = ace(dimor_log, cat.tree)$ace
dimor_ancs_nolog = ace(dimor_nolog, cat.tree)$ace

# estimate ancestral dimorphism by trait in ML
log_data_bytrait = averages$Log_ByTrait_Averages
nolog_data_bytrait = averages$NoLog_ByTrait_Averages

log_dimorps_bytrait = data.frame()
nolog_dimorps_bytrait = data.frame()
for(i in 1:length(log_data_bytrait)){
  # no log
  dimorp_sp = unlist(nolog_data_bytrait[[i]][1]) - unlist(nolog_data_bytrait[[i]][2])
  nolog_dimorps_bytrait = rbind(nolog_dimorps_bytrait, dimorp_sp)
  
  # log
  dimorp_sp = unlist(log_data_bytrait[[i]][1]) - unlist(log_data_bytrait[[i]][2])
  log_dimorps_bytrait = rbind(log_dimorps_bytrait, dimorp_sp)
}

# rename columns and insert species names
# no log
colnames(nolog_dimorps_bytrait) = names(averages$Log_ByTrait_Averages$Allenopithecus_nigroviridis$Machos)
nolog_dimorps_bytrait = cbind(averages$Species, nolog_dimorps_bytrait)
colnames(nolog_dimorps_bytrait)[1] = "Species"

# log
colnames(log_dimorps_bytrait) = names(averages$Log_ByTrait_Averages$Allenopithecus_nigroviridis$Machos)
log_dimorps_bytrait = cbind(averages$Species, log_dimorps_bytrait)
colnames(log_dimorps_bytrait)[1] = "Species"

# ACE
nolog_ancs = sapply(nolog_dimorps_bytrait[,2:40], function(x) ace(x, cat.tree)$ace)
log_ancs = sapply(log_dimorps_bytrait[,2:40], function(x) ace(x, cat.tree)$ace)

nolog_ancs = cbind(seq(119,235,1), nolog_ancs)
colnames(nolog_ancs)[1] = "Species"

log_ancs = cbind(seq(119,235,1), log_ancs)
colnames(log_ancs)[1] = "Species"

# estimate VCV and PCs for ancestral
# read vcv matrices
# read all log vcv matrices
setwd("~/Dropbox/Doc/Output/log_vcv")
temp = list.files(pattern = "*.txt")
log_vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(log_vcv)  = gsub(".txt", replacement= "", temp)

# read all no log vcv matrices
setwd("~/Dropbox/Doc/Output/nolog_vcv")
temp = list.files(pattern = "*.txt")
vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(vcv)  = gsub(".txt", replacement= "", temp)

nolog_all_cov_matrices = PhyloW(cat.tree, vcv)
log_all_cov_matrices = PhyloW(cat.tree, log_vcv)

# PCs of the ancestral species
nolog_pcs = list()
log_pcs = list()
nolog_vals = list()
log_vals = list()
nolog_diags = list()
log_diags = list()
for(i in 119:235){
  # NO LOG
  sp = which(names(nolog_all_cov_matrices) == i)
  covar = nolog_all_cov_matrices[[sp]]
  
  # get 6 PCs and eingenvalues
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
  
  # Diag
  v = diag(as.matrix(covar))
  
  nolog_pcs[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6)
  nolog_vals[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  nolog_diags[[i]] = v
  
  # LOG
  sp = which(names(log_all_cov_matrices) == i)
  covar = log_all_cov_matrices[[sp]]
  
  # get 6 PCs and eingenvalues
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
  
  # Diag
  v = diag(as.matrix(covar))
  
  log_pcs[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6)
  log_vals[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  log_diags[[i]] = v
  
}

# make list and rename it
ancestral_dimorp_pcs = list(dimor_nolog, dimor_log, dimor_ancs_nolog, dimor_ancs_log, nolog_ancs, log_ancs, 
                            nolog_all_cov_matrices[seq(125, 241, 1)], log_all_cov_matrices[seq(125, 241, 1)], 
                            nolog_pcs, log_pcs, nolog_vals, log_vals, nolog_diags, log_diags)

names(ancestral_dimorp_pcs) = c("NoLog_Extant_Dimorphism", "Log_Extant_Dimorphism", "NoLog_Ancestral_Dimorphism",
                                "Log_Ancestral_Dimorphism", "NoLog_ByTrait_Ancestral_Dimorphism", "Log_ByTrait_Ancestral_Dimorphism", 
                                "NoLog_Ancestral_VCV", "Log_Ancestral_VCV", "NoLog_Ancestral_PCs", "Log_Ancestral_PCs",
                                "NoLog_Autovalues", "Log_Autovalues", "NoLog_Diagonal", "Log_Diagonal")

names(ancestral_dimorp_pcs[[1]]) = averages$Species[-index]
names(ancestral_dimorp_pcs[[2]]) = averages$Species[-index]

ancestral_dimorp_pcs$Log_Ancestral_PCs = ancestral_dimorp_pcs$Log_Ancestral_PCs[-seq(1,118,1)]
ancestral_dimorp_pcs$NoLog_Ancestral_PCs = ancestral_dimorp_pcs$NoLog_Ancestral_PCs[-seq(1,118,1)]
names(ancestral_dimorp_pcs$Log_Ancestral_PCs) = names(ancestral_dimorp_pcs$Log_Ancestral_Dimorphism)
names(ancestral_dimorp_pcs$NoLog_Ancestral_PCs) = names(ancestral_dimorp_pcs$NoLog_Ancestral_Dimorphism)


ancestral_dimorp_pcs$Log_Autovalues = ancestral_dimorp_pcs$Log_Autovalues[-seq(1,118,1)]
ancestral_dimorp_pcs$NoLog_Autovalues = ancestral_dimorp_pcs$NoLog_Autovalues[-seq(1,118,1)]
names(ancestral_dimorp_pcs$NoLog_Autovalues) = names(ancestral_dimorp_pcs$NoLog_Ancestral_Dimorphism)
names(ancestral_dimorp_pcs$Log_Autovalues) = names(ancestral_dimorp_pcs$NoLog_Ancestral_Dimorphism)

ancestral_dimorp_pcs$Log_Diagonal = ancestral_dimorp_pcs$Log_Diagonal[-seq(1,118,1)]
ancestral_dimorp_pcs$NoLog_Diagonal = ancestral_dimorp_pcs$NoLog_Diagonal[-seq(1,118,1)]
names(ancestral_dimorp_pcs$Log_Diagonal) = names(ancestral_dimorp_pcs$NoLog_Ancestral_Dimorphism)
names(ancestral_dimorp_pcs$NoLog_Diagonal) = names(ancestral_dimorp_pcs$NoLog_Ancestral_Dimorphism)

for(i in 1:117){
  names(ancestral_dimorp_pcs$Log_Ancestral_PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  names(ancestral_dimorp_pcs$NoLog_Ancestral_PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  names(ancestral_dimorp_pcs$Log_Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6")
  names(ancestral_dimorp_pcs$NoLog_Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6")
}

# save final list of results
saveRDS(ancestral_dimorp_pcs, file = "~/Dropbox/Novo_Output/Averages_PCS_Ancestral_Species.RData")
