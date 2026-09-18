# Average and PCs of extant species
# log and no-log averages and PCs

##################################################################################
### Catarrhini
##################################################################################

# load packages
if(!require(stringr)){install.packages("stringr"); library(stringr)}

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

# read all no log vcv matrices
setwd("~/Dropbox/Doc/Output/nolog_vcv")
temp = list.files(pattern = "*.txt")

# remove species which are not in the phylogeny
to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), temp)
temp = temp[-index]

vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(vcv)  = gsub(".txt", replacement= "", temp)

# read all log vcv matrices
setwd("~/Dropbox/Doc/Output/log_vcv")
temp = list.files(pattern = "*.txt")

# remove species which are not in the phylogeny
to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), temp)
temp = temp[-index]

log_vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(log_vcv)  = gsub(".txt", replacement= "", temp)

# create lists
averages_extant_log = list()
averages_extant_nolog = list()
bytrait_extant_log = list()
bytrait_extant_nolog = list()
pcs_extant_log = list()
pcs_extant_nolog = list()
vals_extant_log = list()
vals_extant_nolog = list()
diags = list()
log_diags = list()
for(i in 1:length(names(vcv))){
  # choose species i
  genus = str_split_1(names(vcv)[[i]], "_")[1]
  sp = str_split_1(names(vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES == sp), ]
  
  # separate M and F
  sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]
  
  # NO LOG
  # geometric means of M and F
  avg_m = apply(sub_sexes_m[,49:87], 2, geomean)
  avg_f = apply(sub_sexes_f[,49:87], 2, geomean)
  
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
  
  # Diag
  v = diag(as.matrix(covar))

  # insert in list of species i
  averages_extant_nolog[[i]] = list(gen_avg_m, gen_avg_f)
  bytrait_extant_nolog[[i]] = list(avg_m, avg_f)
  pcs_extant_nolog[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6)
  vals_extant_nolog[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  diags[[i]] = v
  
  # LOG
  # geometric means of M and F
  avg_m = apply(log(sub_sexes_m[,49:87] * 10), 2, geomean)
  avg_f = apply(log(sub_sexes_f[,49:87] * 10), 2, geomean)
  
  # general geometric means of M and F
  gen_avg_m = geomean(avg_m)
  gen_avg_f = geomean(avg_f)
  
  # PCs ans eigenvalues
  covar = log_vcv[[i]]
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
  v_log = diag(as.matrix(covar))
  
  # insert in list of species i
  averages_extant_log[[i]] = list(gen_avg_m, gen_avg_f)
  bytrait_extant_log[[i]] = list(avg_m, avg_f)
  pcs_extant_log[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6)
  vals_extant_log[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  log_diags[[i]] = v_log
}

# get names of species
especies = names(vcv)

# create and naming final list
extant_averages_pcs = list(especies, averages_extant_log, averages_extant_nolog, bytrait_extant_log, 
                    bytrait_extant_nolog, pcs_extant_log, pcs_extant_nolog, vals_extant_nolog,
                    vals_extant_log, diags, log_diags)

names(extant_averages_pcs) = c("Species", "Log_Averages", "NoLog_Averages", "Log_ByTrait_Averages",
                        "NoLog_ByTrait_Averages", "Log_PCs", "NoLog_PCs", "NoLog_Autovalues", 
                        "Log_Autovalues", "NoLog_Diagonal", "Log_Diagonal")

names(extant_averages_pcs$Log_Averages) = especies
names(extant_averages_pcs$NoLog_Averages) = especies
names(extant_averages_pcs$Log_ByTrait_Averages) = especies
names(extant_averages_pcs$NoLog_ByTrait_Averages) = especies
names(extant_averages_pcs$Log_PCs) = especies
names(extant_averages_pcs$NoLog_PCs) = especies
names(extant_averages_pcs$NoLog_Autovalues) = especies
names(extant_averages_pcs$Log_Autovalues) = especies
names(extant_averages_pcs$NoLog_Diagonal) = especies
names(extant_averages_pcs$Log_Diagonal) = especies

for(i in 1:118){
  names(extant_averages_pcs$Log_Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs$NoLog_Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs$Log_ByTrait_Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs$NoLog_ByTrait_Averages[[i]]) = c("Machos", "Fêmeas")
}

for(i in 1:118){
  names(extant_averages_pcs$Log_PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  names(extant_averages_pcs$NoLog_PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  names(extant_averages_pcs$NoLog_Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda", "Lambda5", "Lambda6")
  names(extant_averages_pcs$Log_Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6")
}

# save the final list in Output folder
#saveRDS(extant_averages_pcs, file = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Averages_PCS_Extant_Species.RData")


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
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini.csv", dec = ",", sep = ",")

# new column with species names
species = vector()
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUB.[i])){
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], msrs$SUB.[i], sep="_")
  }
}

msrs = cbind(species, msrs)
species = unique(species)

index = which(apply(msrs[,24:62], 1, anyNA) == TRUE)
msrs = msrs[-index,]

# remove uncertain sex
msrs$SEX4.[msrs$SEX4. == ""] = NA
msrs$SEX4.[msrs$SEX4. == " "] = NA

msrs = msrs[complete.cases(msrs$SEX4.), ]

# read all no log vcv matrices
setwd("~/Dropbox/Doc/Output/p_nolog_vcv")
temp = list.files(pattern = "*.txt")
vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(vcv)  = gsub(".txt", replacement= "", temp)

# read all log vcv matrices
setwd("~/Dropbox/Doc/Output/p_log_vcv")
temp = list.files(pattern = "*.txt")
log_vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(log_vcv)  = gsub(".txt", replacement= "", temp)

# create lists
averages_extant_log = list()
averages_extant_nolog = list()
bytrait_extant_log = list()
bytrait_extant_nolog = list()
pcs_extant_log = list()
pcs_extant_nolog = list()
vals_extant_log = list()
vals_extant_nolog = list()
diags = list()
log_diags = list()
for(i in 1:length(names(vcv))){
  # choose species i
  genus = str_split_1(names(vcv)[[i]], "_")[1]
  sp = str_split_1(names(vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS. == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES. == sp), ]
  
  # separate M and F
  sub_sexes_m = species_subset[which(species_subset$SEX4. == "M"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX4. == "F"), ]
  
  # NO LOG
  # geometric means of M and F
  avg_m = apply(sub_sexes_m[,24:62], 2, geomean)
  avg_f = apply(sub_sexes_f[,24:62], 2, geomean)
  
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
  
  # Diag
  v = diag(as.matrix(covar))
  
  # insert in list of species i
  averages_extant_nolog[[i]] = list(gen_avg_m, gen_avg_f)
  bytrait_extant_nolog[[i]] = list(avg_m, avg_f)
  pcs_extant_nolog[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6)
  vals_extant_nolog[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  diags[[i]] = v
  
  # LOG
  # geometric means of M and F
  avg_m = apply(log(sub_sexes_m[,24:62] * 10), 2, geomean)
  avg_f = apply(log(sub_sexes_f[,24:62] * 10), 2, geomean)
  
  # general geometric means of M and F
  gen_avg_m = geomean(avg_m)
  gen_avg_f = geomean(avg_f)
  
  # PCs ans eigenvalues
  covar = log_vcv[[i]]
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
  v_log = diag(as.matrix(covar))
  
  # insert in list of species i
  averages_extant_log[[i]] = list(gen_avg_m, gen_avg_f)
  bytrait_extant_log[[i]] = list(avg_m, avg_f)
  pcs_extant_log[[i]] = list(pc_1, pc_2, pc_3, pc_4, pc_5, pc_6)
  vals_extant_log[[i]] = list(vals_1, vals_2, vals_3, vals_4, vals_5, vals_6)
  log_diags[[i]] = v_log
}

# get names of species
especies = names(vcv)

# create and naming final list
extant_averages_pcs = list(especies, averages_extant_log, averages_extant_nolog, bytrait_extant_log, 
                           bytrait_extant_nolog, pcs_extant_log, pcs_extant_nolog, vals_extant_nolog,
                           vals_extant_log, diags, log_diags)

names(extant_averages_pcs) = c("Species", "Log_Averages", "NoLog_Averages", "Log_ByTrait_Averages",
                               "NoLog_ByTrait_Averages", "Log_PCs", "NoLog_PCs", "NoLog_Autovalues", 
                               "Log_Autovalues", "NoLog_Diagonal", "Log_Diagonal")

names(extant_averages_pcs$Log_Averages) = especies
names(extant_averages_pcs$NoLog_Averages) = especies
names(extant_averages_pcs$Log_ByTrait_Averages) = especies
names(extant_averages_pcs$NoLog_ByTrait_Averages) = especies
names(extant_averages_pcs$Log_PCs) = especies
names(extant_averages_pcs$NoLog_PCs) = especies
names(extant_averages_pcs$NoLog_Autovalues) = especies
names(extant_averages_pcs$Log_Autovalues) = especies
names(extant_averages_pcs$NoLog_Diagonal) = especies
names(extant_averages_pcs$Log_Diagonal) = especies

for(i in 1:72){
  names(extant_averages_pcs$Log_Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs$NoLog_Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs$Log_ByTrait_Averages[[i]]) = c("Machos", "Fêmeas")
  names(extant_averages_pcs$NoLog_ByTrait_Averages[[i]]) = c("Machos", "Fêmeas")
}

for(i in 1:72){
  names(extant_averages_pcs$Log_PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  names(extant_averages_pcs$NoLog_PCs[[i]]) = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  names(extant_averages_pcs$NoLog_Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda", "Lambda5", "Lambda6")
  names(extant_averages_pcs$Log_Autovalues[[i]]) = c("Lambda1", "Lambda2", "Lambda3", "Lambda4", "Lambda5", "Lambda6")
}

# save the final list in Output folder
#saveRDS(extant_averages_pcs, file = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/p_Averages_PCS_Extant_Species.RData")

