# set WD and functions
setwd("~/Dropbox/Doc/Data/wos_mating_systems/")

# geometric mean and vector correlation functions
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# Read and unique species
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini.csv", dec = ",", sep = ",")

msrs$GENUS. =  gsub(" ", "", msrs$GENUS.)

msrs$SPECIES.[which(msrs$SPECIES. == "flaviceps?")] = "flaviceps"
msrs$SPECIES.[which(msrs$SPECIES. == "xanthosternos")] = "xanthosternus"

msrs$SEX.[msrs$SEX. == "sem"] = NA
msrs$SEX.[msrs$SEX. == "Sem"] = NA
msrs$SEX.[msrs$SEX. == "sem "] = NA
msrs$SEX.[msrs$SEX. == "sem?"] = NA
msrs$SEX.[msrs$SEX. == "F?"] = NA
msrs$SEX.[msrs$SEX. == ""] = NA
msrs$SEX.[msrs$SEX. == "CADE?"] = NA
msrs$SEX.[msrs$SEX. == "?"] = NA

msrs = msrs[complete.cases(msrs$SEX.), ]

msrs$SUB.[which(msrs$SUB. == "")] = "NA"

species = vector()
# new column with species names
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUB.[i])){
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], msrs$SUB.[i], sep="_")
  }
}

msrs = cbind(species, msrs)
species = unique(species)




#read V/CV matrices
#setwd("~/Dropbox/Doc/Data/matriesvcv/")

# read all V/CV matrices
#temp = list.files(pattern="*.txt")
#vcv = lapply(temp, read.table)
#names(vcv)  = gsub(".txt", replacement= "", temp)
#names(vcv) = gsub("VCV2", "", names(vcv))

# loop to get Dimorphisms, PCs and vector correlations
#dimorp_pcs = data.frame()
#for(i in 1:length(species)){
#  #for each genus in V/CV matrices
#  genero = sub("\\_.*", "", species[i])
#  index = grepl(genero, names(vcv), fixed = TRUE)
#  index = as.integer(as.logical(index))
#  covar = vcv[[which(index == 1)]]
  
  #covar = vcv[["Hylobates"]]
  #x = setNames(pc_4, names(msrs[,50:88]))
  
  #x[x > 0]
  #x[x < 0]
  
  # PCs
  #pc_1 = eigen(covar)$vectors[,1]
  #pc_2 = eigen(covar)$vectors[,2]
  #pc_3 = eigen(covar)$vectors[,3]
  #pc_4 = eigen(covar)$vectors[,4]
  
  # correlations between isometric vector and PCs 1 and 2
  #cor_iso_1 = abs(corVector(pc_1, rep(0.160128154, length(pc_1))))
  #cor_iso_2 = abs(corVector(pc_2, rep(0.160128154, length(pc_2))))
  #cor_iso_3 = abs(corVector(pc_3, rep(0.160128154, length(pc_3))))
  #cor_iso_4 = abs(corVector(pc_4, rep(0.160128154, length(pc_4))))
  
  # choose species
  sub_species = msrs[which(msrs$species == species[1]), ]
  
  # separate M and F
  sub_sexes_m = sub_species[which(sub_species$SEX. == "M"), ]
  sub_sexes_f = sub_species[which(sub_species$SEX. == "F"), ]
  
  # sexual dimorphism with geometric means
  avg_m = apply(sub_sexes_m[,24:62], 2, geomean)
  avg_f = apply(sub_sexes_f[,24:62], 2, geomean)
  
  # geom average for males and females
  gen_avg_m = geomean(avg_m)
  gen_avg_f = geomean(avg_f)
  
  # Normalized Sexual Dimorphisms
  ind_dimorp = (avg_m - avg_f) / gen_avg_m
  #ind_dimorp = avg_m - avg_f
  
  # correlations between dimorphisms and PCs
  #c_1 = abs(corVector(ind_dimorp, pc_1))
  #c_2 = abs(corVector(ind_dimorp, pc_2))
  #c_3 = abs(corVector(ind_dimorp, pc_3))
  #c_4 = abs(corVector(ind_dimorp, pc_4))
  
  # norms of dimorphism vectors
  norms = norma(ind_dimorp)
  
  # final data frame
  cors = c(c_1, c_2, c_3, c_4, cor_iso_1, cor_iso_2, cor_iso_3, cor_iso_4, norms)
  dimorp_pcs = rbind(dimorp_pcs, cors)
  
}

# insert names of species and rearrange data frame
colnames(dimorp_pcs) = c("dimorp_PC1", "dimorp_PC2", "dimorp_PC3", "dimorp_PC4", "cor_iso_pc1", "cor_iso_pc2", "cor_iso_pc3", "cor_iso_pc4", "Dimorp_Norms")
dimorp_pcs = cbind(dimorp_pcs, species)
dimorp_pcs = dimorp_pcs[, c("species", "dimorp_PC1", "dimorp_PC2", "dimorp_PC3", "dimorp_PC4", "cor_iso_pc1", "cor_iso_pc2", "cor_iso_pc3", "cor_iso_pc4", "Dimorp_Norms")]

# read mating systems table
matings = read.table("~/Dropbox/Doc/Data/wos_mating_systems/primates_mating_systems.csv",
                     sep = ",", header = TRUE)

# get only the species of "interest"
species = vector()
for(i in 1:nrow(matings)){
  if(is.na(matings$SUBSPECIES[i])){
    species[i] = paste(matings$GENUS[i], matings$SPECIES[i], sep="_")}
  else{
    species[i] = paste(matings$GENUS[i], matings$SPECIES[i], matings$SUBSPECIES[i], sep="_")
  }
}

# prepare final data frame
final = data.frame()
for(i in 1:length(species)){
  index = which(species[i] == dimorp_pcs$species)
  final = rbind(final, dimorp_pcs[index,])
  
}

# remove _ and insert spaces.
final$species = sub("_", " ", final$species)

# save table
#write.table(final, file = "~/Dropbox/Doc/Output/Pcs_Dimorp_Norma.txt", row.names = FALSE, 
#            dec = ".", sep = '\t', quote = FALSE)

# plots and correlations (important result)
par(mfrow = c(2, 2))
plot(final$dimorp_PC1, final$Dimorp_Norms, xlab = "Correlação Vetor de Dimorfismo e PC1", ylab = "Norma do Vetor de Dimorfismo")
plot(final$dimorp_PC2, final$Dimorp_Norms, xlab = "Correlação Vetor de Dimorfismo e PC2", ylab = "Norma do Vetor de Dimorfismo")
plot(final$dimorp_PC3, final$Dimorp_Norms, xlab = "Correlação Vetor de Dimorfismo e PC3", ylab = "Norma do Vetor de Dimorfismo")
plot(final$dimorp_PC4, final$Dimorp_Norms, xlab = "Correlação Vetor de Dimorfismo e PC4", ylab = "Norma do Vetor de Dimorfismo")

corVector(final$dimorp_PC1, final$Dimorp_Norms)
corVector(final$dimorp_PC2, final$Dimorp_Norms)
corVector(final$dimorp_PC3, final$Dimorp_Norms)
corVector(final$dimorp_PC4, final$Dimorp_Norms)
