#Betas and DeltaZs
# log and no log

#setwd

#load packages and phylogeny
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(ape)){install.packages("ape"); library(ape)}

filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)

# read species names
setwd("~/Dropbox/Doc/Output/nolog_vcv")
temp = list.files(pattern = "*.txt")
names(temp)  = gsub(".txt", replacement= "", temp)

# remove species which are not in the phylogeny
to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), names(temp))
temp = temp[-index]
cat.tree = keep.tip(tree, names(temp))

# load extant and ancestral average dimorphism
extant = readRDS("~/Dropbox/Novo_Output/Averages_PCS_Extant_Species.RData")
ancestral = readRDS("~/Dropbox/Novo_Output/Averages_PCS_Ancestral_Species.RData")

# get the dimorphism for all species and nodes in phylogeny
log_dimor_ancestral = ancestral$Log_ByTrait_Ancestral_Dimorphism
nolog_dimor_ancestral = ancestral$Log_ByTrait_Ancestral_Dimorphism

nolog_dimor_extant = list()
log_dimor_extant = list()
for(i in 1:length(extant$Species)){
  nolog_dimor_extant[[i]] = unlist(extant$NoLog_ByTrait_Averages[[i]][1]) - unlist(extant$NoLog_ByTrait_Averages[[i]][2])
  log_dimor_extant[[i]] = unlist(extant$Log_ByTrait_Averages[[i]][1]) - unlist(extant$Log_ByTrait_Averages[[i]][2])
}

# name the dimorphisms
names(nolog_dimor_extant) = extant$Species
names(log_dimor_extant) = extant$Species

# NO LOG
# give node names for deltaZ calculation
node_names = c(cat.tree$tip.label, 119:235)
node_numbers = c(1:118, 120:235)
nolog_delta_Zs = vector("list", length(node_numbers))

lista_medias = list(c(nolog_dimor_extant, apply(nolog_dimor_ancestral[,2:40], 1, as.list)))

# calculate Delta_Z
for(i in 1:length(node_numbers)){
  node = node_numbers[i]
  current_node = node_names[node]
  node_mean = lista_medias[[1]][[current_node]]
  ancestral_node = node_names[cat.tree$edge[cat.tree$edge[,2] == node, 1]]
  ancestral_mean = lista_medias[[1]][[ancestral_node]] 
  nolog_delta_Zs[[i]] = unlist(node_mean) - unlist(ancestral_mean)
}

names(nolog_delta_Zs) = node_names[-119]
nolog_lista_deltaZs = ldply(nolog_delta_Zs)

# log
# give node names for deltaZ calculation
node_names = c(cat.tree$tip.label, 119:235)
node_numbers = c(1:118, 120:235)
log_delta_Zs = vector("list", length(node_numbers))

lista_medias = list(c(log_dimor_extant, apply(log_dimor_ancestral[,2:40], 1, as.list)))

# calculate Delta_Z
for(i in 1:length(node_numbers)){
  node = node_numbers[i]
  current_node = node_names[node]
  node_mean = lista_medias[[1]][[current_node]]
  ancestral_node = node_names[cat.tree$edge[cat.tree$edge[,2] == node, 1]]
  ancestral_mean = lista_medias[[1]][[ancestral_node]] 
  log_delta_Zs[[i]] = unlist(node_mean) - unlist(ancestral_mean)
}

names(log_delta_Zs) = node_names[-119]
log_lista_deltaZs = ldply(log_delta_Zs)

# Beta calculation
# B = P^-1 deltaz

# read all no log vcv matrices
setwd("~/Dropbox/Doc/Output/nolog_vcv")
temp = list.files(pattern = "*.txt")

to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), temp)
temp = temp[-index]

vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(vcv)  = gsub(".txt", replacement= "", temp)

# read all log vcv matrices
setwd("~/Dropbox/Doc/Output/log_vcv")
temp = list.files(pattern = "*.txt")

to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), temp)
temp = temp[-index]

log_vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(log_vcv)  = gsub(".txt", replacement= "", temp)

# Matrix inversion of P with ExtendMatrix of evolqg
i_vcv = list()
i_log_vcv = list()
for(i in 1:length(vcv)){
  i_vcv[[i]] = ExtendMatrix(vcv[[i]])$ExtMat
  i_log_vcv[[i]] = ExtendMatrix(log_vcv[[i]])$ExtMat
}

names(i_vcv) = names(vcv)
names(i_log_vcv) = names(log_vcv)

# Betas
nolog_extant_betas = list()
log_extant_betas = list()
for(i in 1:length(vcv)){
  # extant
  nome_sp = names(vcv)[i]
  
  # no log
  deltaz = nolog_lista_deltaZs[nolog_lista_deltaZs$.id == nome_sp,][2:40] 
  matrizp = i_vcv[[nome_sp]]
  nolog_extant_betas[[i]] = solve(matrizp, deltaz)
  names(nolog_extant_betas)[[i]] = nome_sp
      
  # log
  deltaz = log_lista_deltaZs[nolog_lista_deltaZs$.id == nome_sp,][2:40]
  matrizp = i_log_vcv[[nome_sp]]
  log_extant_betas[[i]] = solve(matrizp, deltaz)
  names(log_extant_betas)[[i]] = nome_sp
}

nolog_ancestral_betas = list()
log_ancestral_betas = list()
for(i in 1:dim(ancestral$NoLog_ByTrait_Ancestral_Dimorphism)[1]){
  # ancestral
  nome_sp = names(ancestral$NoLog_Ancestral_Dimorphism)[i]
  
  # no log
  deltaz = nolog_lista_deltaZs[nolog_lista_deltaZs$.id == nome_sp,][2:40]
  matrizp = ExtendMatrix(ancestral$NoLog_Ancestral_VCV[[nome_sp]])$ExtMat 
  nolog_ancestral_betas[[i]] = solve(matrizp, deltaz)
  names(nolog_ancestral_betas)[[i]] = nome_sp
  
  # log
  deltaz = log_lista_deltaZs[nolog_lista_deltaZs$.id == nome_sp,][2:40]
  matrizp = ExtendMatrix(ancestral$Log_Ancestral_VCV[[nome_sp]])$ExtMat
  log_ancestral_betas[[i]] = solve(matrizp, deltaz)
  names(log_ancestral_betas)[[i]] = nome_sp
}


# form a single list
betas_deltaz = list(nolog_lista_deltaZs, log_lista_deltaZs, nolog_extant_betas, log_extant_betas,
                    nolog_ancestral_betas, log_ancestral_betas)

# rename the list components
names(betas_deltaz) = c("NoLog_DeltaZs", "Log_DeltaZs", "NoLog_Extant_Betas", "Log_Extant_Betas",
                        "NoLog_Ancestral_Betas", "Log_Ancestral_Betas")

# save the final list
saveRDS(betas_deltaz, file = "~/Dropbox/Novo_Output/Betas_DeltaZ.RData")
