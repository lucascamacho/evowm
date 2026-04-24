# Average and PCs of extant species
##################################################################################
### Catarrhini
##################################################################################

setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")

# load packages
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(ape)){install.packages("ape"); library(ape)}

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# Read and unique species
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")

# check names and remove doubts
msrs$SEX[which(msrs$SEX == "?female")] = "female"
msrs$SEX[which(msrs$SEX == "?male")] = "male"

# remove uncertain sex
msrs$SEX[msrs$SEX == "0"] = NA
msrs$SEX[msrs$SEX == ""] = NA
msrs$SEX[msrs$SEX == "sexo"] = NA

msrs = msrs[complete.cases(msrs$SEX), ]

# read all vcv matrices
#setwd("~/Dropbox/Doc/Data/p_vcv_gabriel/catarrhini")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
names(vcv)  = gsub(".csv", replacement= "", temp)

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
species = names(vcv)
tree = drop.tip(tree, setdiff(tree$tip.label, species))
species <- species[match(tree$tip.label, species)]

# create lists
averages = list()
bytrait_averages = list()
pcs = list()
vals = list()
diags = list()

n_cata = length(species[1:49])
for(i in 1:n_cata){
  # choose species i
  genus = str_split_1(species[i], "_")[1]
  sp = str_split_1(species[i], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES == sp), ]

    # separate M and F
  sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]
  
  # geometric means of M and F
  avg_m = apply(log(sub_sexes_m[, c(49:67, 69:87)]), 2, geomean)
  avg_f = apply(log(sub_sexes_f[, c(49:67, 69:87)]), 2, geomean)

  
  # general geometric means of M and F
  gen_avg_m = geomean(avg_m)
  gen_avg_f = geomean(avg_f)
  
  # PCs and eigenvalues
  covar = vcv[[i]]
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
  
  # Diag
  v = diag(as.matrix(covar))
  
  # insert in list of species i
  averages[[i]] = list(gen_avg_m, gen_avg_f)
  bytrait_averages[[i]] = list(avg_m, avg_f)
  pcs[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6, pc_7, pc_8)
  vals[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  diags[[i]] = v
}
# get names of species
especies = species[1:n_cata]

# create and naming final list
extant_averages_pcs_c = list(especies, averages, bytrait_averages, pcs, vals, diags)

names(extant_averages_pcs_c) = c("Species", "Averages", "ByTrait_Averages", "PCs", "Autovalues", 
                               "Diagonal")

names(extant_averages_pcs_c$Averages) = especies
names(extant_averages_pcs_c$ByTrait_Averages) = especies
names(extant_averages_pcs_c$PCs) = especies
names(extant_averages_pcs_c$Autovalues) = especies
names(extant_averages_pcs_c$Diagonal) = especies

for(i in 1:length(especies)){
  names(extant_averages_pcs_c$Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs_c$ByTrait_Averages[[i]]) = c("Machos", "Fêmeas")
}

for(i in 1:length(especies)){
  names(extant_averages_pcs_c$PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  names(extant_averages_pcs_c$Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda", "Lambda5", "Lambda6")
}

##################################################################################
### Platyrrhini
##################################################################################

# load packages
if(!require(stringr)){install.packages("stringr"); library(stringr)}

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# Read and unique species
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini.csv", dec = ".", sep = ",")

index = which(apply(msrs[, 23:61], 1, function(x) any(is.na(x) | (x == "" & !is.numeric(x)))))
msrs = msrs[-index,]

# remove uncertain sex
msrs$SEX4.[msrs$SEX4. == ""] = NA
msrs$SEX4.[msrs$SEX4. == " "] = NA

msrs = msrs[complete.cases(msrs$SEX4.), ]

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
species = names(vcv)
tree = drop.tip(tree, setdiff(tree$tip.label, species))
species <- species[match(tree$tip.label, species)]

# create lists
averages= list()
bytrait_averages = list()
pcs = list()
vals = list()
diags = list()

n_platy = length(species[50:73])
species = species[50:73]
for(i in 1:n_platy){
  # choose species i
  genus = str_split_1(species[i], "_")[1]
  sp = str_split_1(species[i], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS. == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES. == sp), ]

  # separate M and F
  sub_sexes_m = species_subset[which(species_subset$SEX4. == "M"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX4. == "F"), ]
  
  # geometric means of M and F
  avg_m = apply(log(sub_sexes_m[,c(23:41,  43:61)]), 2, geomean)
  avg_f = apply(log(sub_sexes_f[,c(23:41,  43:61)]), 2, geomean)
  
  # general geometric means of M and F
  gen_avg_m = geomean(avg_m)
  gen_avg_f = geomean(avg_f)
  
  # PCs and eigenvalues
  covar = vcv[[i]]
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
  
  # Diag
  v = diag(as.matrix(covar))
  
  # insert in list of species i
  extant_averages_pcs_c$Averages[[49+i]] = list(gen_avg_m, gen_avg_f)
  extant_averages_pcs_c$ByTrait_Averages[[49+i]] = list(avg_m, avg_f)
  extant_averages_pcs_c$PCs[[49+i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6, pc_7, pc_8)
  extant_averages_pcs_c$Autovalues[[49+i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  extant_averages_pcs_c$Diagonal[[49+i]] = v
}

# get names of species
extant_averages_pcs_c$Species = append(extant_averages_pcs_c$Species, species)

#
names(extant_averages_pcs_c$Averages)[50:73] = species
names(extant_averages_pcs_c$ByTrait_Averages)[50:73] = species
names(extant_averages_pcs_c$PCs)[50:73] = species
names(extant_averages_pcs_c$Autovalues)[50:73] = species
names(extant_averages_pcs_c$Diagonal)[50:73] = species

#
for(i in 50:73){
  names(extant_averages_pcs_c$Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs_c$ByTrait_Averages[[i]]) = c("Machos", "Fêmeas")
}

#
for(i in 50:73){
  names(extant_averages_pcs_c$PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")
  names(extant_averages_pcs_c$Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda", "Lambda5", "Lambda6")
}

saveRDS(extant_averages_pcs_c, "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")
