# set WD, packages and functions
setwd("~/Dropbox/Doc/Data/wos_mating_systems/")

library(ggplot2)
library(dplyr)
library(vegan)
require(plyr)

Classify = function(vector){
  drn = factor(vector)
  return(as.numeric(drn))
}

geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# read and unique species
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")
species = vector()

msrs$SUBSPECIES[which(msrs$SUBSPECIES == "0")] = NA
msrs$SUBSPECIES[which(msrs$SUBSPECIES == "")] = NA


msrs$SEX[which(msrs$SEX == "?female")] = "female"
msrs$SEX[which(msrs$SEX == "?male")] = "male"

msrs$SEX[msrs$SEX == "0"] = NA
msrs$SEX[msrs$SEX == ""] = NA
msrs$SEX[msrs$SEX == "sexo"] = NA

msrs = msrs[complete.cases(msrs$SEX), ]

# merge info of species in sigle cell
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUBSPECIES[i])){
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], msrs$SUBSPECIES[i], sep="_")
  }
}

msrs = cbind(species, msrs)
species = unique(species)

# get averages, index of dimorphism
avg_male = vector()
avg_female = vector()
ind_dimorp = vector()

for(i in 1:length(species)){
  # choose species
  sub_species = msrs[which(msrs$species == species[i]), ]
  
  # separate in M and F
  sub_sexes_m = sub_species[which(sub_species$SEX == "male"), ]
  sub_sexes_f = sub_species[which(sub_species$SEX == "female"), ]
  
  # Averages of M and F
  avg_m = apply(sub_sexes_m[,50:88], 2, geomean)
  avg_f = apply(sub_sexes_f[,50:88], 2, geomean)
  
  avg_male[i] = geomean(avg_m)
  avg_female[i] = geomean(avg_f)
  
  # indice of dimorphism
  ind_dimorp[i] = geomean(avg_m) / geomean(avg_f)

}

data = data.frame(species, avg_male, avg_female, ind_dimorp)
data = na.omit(data)

# Mating data
dados = read.csv("primates_mating_systems.csv", header = TRUE)
dados$PROP_MALES_FEMALES =  as.numeric(gsub(",", ".", gsub("\\.", "", dados$PROP_MALES_FEMALES)))

# classify mating data to numbers
dados[5:6] = apply(dados[5:6], 2, Classify)
dados[8:13] = apply(dados[8:13], 2, Classify)

# merge species names in a single cell
for(i in 1:nrow(dados)){
  if(is.na(dados$SUBSPECIES[i])){
    dados$especie[i] = paste(dados$GENUS[i], dados$SPECIES[i], sep="_")}
  else{
    dados$especie[i] = paste(dados$GENUS[i], dados$SPECIES[i], dados$SUBSPECIES[i], sep="_")
  }
}

# MDS
d = dist(dados[ ,5:9])
fit = cmdscale(d, eig = TRUE, k = 2)

x = fit$points[,1]
y = fit$points[,2]

# get only the "interested" species
final = data.frame()
for(i in 1:nrow(dados)){
  index = which(dados$especie[i] == data$species)
  final = rbind(final, data[index,])
  
}

# get all data and rename data.frame
paraplot = cbind(final, x, y, dados$SOCIAL_ORGANIZATION, dados$MATING_SYSTEM, dados$PROP_MALES_FEMALES, 
             dados$DOMINANCE, dados$AGGRESSION)

colnames(paraplot) = c("species", "avg_male", "avg_female", "ind_dimorphism", "mds1", "mds2", 
                       "SOCIAL_ORGANIZATION", "MATING_SYSTEM", "PROP_MALES_FEMALES", "DOMINANCE", "AGGRESSION")
row.names(paraplot) = NULL

paraplot$species = sub("_", " ", paraplot$species)

# saving table with "," decimal for Mesquite only
#paraplot[,2:11] = numcolwise(prettyNum)(paraplot[,2:11], dec = ",")

#write.table(paraplot, file = "~/Dropbox/Doc/Output/Dimorp_MDS_OrgVars.txt.txt", row.names = FALSE, dec = ".", 
#            sep = '\t', quote = FALSE)

