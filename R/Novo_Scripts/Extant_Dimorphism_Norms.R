#Dimorphism vector norms
# log and no-log dimorphisms
# set WD

# load packages
if(!require(stringr)){install.packages("stringr"); library(stringr)}

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# vector correlation vector
prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

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

# read all species names and remove species out of the phylogeny
setwd("~/Dropbox/Doc/Output/nolog_vcv")
vcv = list.files(pattern = "*.txt")

# remove species which are not in the phylogeny
to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), temp)
vcv = vcv[-index]
names(vcv)  = gsub(".txt", replacement= "", vcv)

log_norma = list()
nolog_norma = list()
log_ind_dimor = list()
nolog_ind_dimor = list()
for(i in 1:length(names(vcv))){
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
  
  # Normalized Sexual Dimorphisms
  ind_dimorp = (avg_m - avg_f) / gen_avg_m
  
  # Norm of dimorphism vectors
  nor = norma(ind_dimorp)
  
  # insert results in list
  nolog_ind_dimor[[i]] = ind_dimorp
  nolog_norma[[i]] = nor
  
  # LOG
  # geometric means of M and F
  avg_m = apply(log(sub_sexes_m[,49:87] * 10), 2, geomean)
  avg_f = apply(log(sub_sexes_f[,49:87] * 10), 2, geomean)
  
  # general geometric means of M and F
  gen_avg_m = geomean(avg_m)
  gen_avg_f = geomean(avg_f)
  
  # Normalized Sexual Dimorphisms
  ind_dimorp = (avg_m - avg_f) / gen_avg_m
  
  # Norm of dimorphism vectors
  nor = norma(ind_dimorp)
  
  # insert results in a list
  log_ind_dimor[[i]] = ind_dimorp
  log_norma[[i]] = nor
  
}

# create and rename final list
normas_dimorfismos = list(log_norma, nolog_norma, log_ind_dimor, nolog_ind_dimor)

names(normas_dimorfismos) = c("Log_Norm_Dimorphism", "NoLog_Norm_Dimorphism", "Log_Normalized_Dimorphism",
                              "NoLog_Normalized_Dimorphism")

names(normas_dimorfismos$Log_Norm_Dimorphism) = names(vcv)
names(normas_dimorfismos$NoLog_Norm_Dimorphism) = names(vcv)
names(normas_dimorfismos$Log_Normalized_Dimorphism) = names(vcv)
names(normas_dimorfismos$NoLog_Normalized_Dimorphism) = names(vcv)

# save final list
saveRDS(normas_dimorfismos, file = "~/Dropbox/Novo_Output/Extant_Dimorphism_Norms.RData")

