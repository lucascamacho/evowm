# set WD and functions
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

# load packages
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(openxlsx)){install.packages("openxlsx"); library(openxlsx)}

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

# Read and unique species
msrs_catarrhini = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")

# Check catarrhini mrsr
# check names and remove doubts
msrs_catarrhini$SEX[which(msrs_catarrhini$SEX == "?female")] = "female"
msrs_catarrhini$SEX[which(msrs_catarrhini$SEX == "?male")] = "male"

# remove uncertain sex
msrs_catarrhini$SEX[msrs_catarrhini$SEX == "0"] = NA
msrs_catarrhini$SEX[msrs_catarrhini$SEX == ""] = NA
msrs_catarrhini$SEX[msrs_catarrhini$SEX == "sexo"] = NA

msrs_catarrhini = msrs_catarrhini[complete.cases(msrs_catarrhini$SEX), ]

#
msrs_platyrrhini = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini.csv", dec = ".", sep = ",")

# Check platyrrhini msrs
index = which(apply(msrs_platyrrhini[, 23:61], 1, function(x) any(is.na(x) | (x == "" & !is.numeric(x)))))
msrs_platyrrhini = msrs_platyrrhini[-index,]

# remove uncertain sex
msrs_platyrrhini$SEX4.[msrs_platyrrhini$SEX4. == ""] = NA
msrs_platyrrhini$SEX4.[msrs_platyrrhini$SEX4. == " "] = NA

msrs_platyrrhini = msrs_platyrrhini[complete.cases(msrs_platyrrhini$SEX4.), ]

#####################################################################################################
# VCV 
# get species names
setwd("~/Dropbox/Doc/Data/p_vcv_gabriel/catarrhini")
temp = list.files(pattern = "*.csv")
species  = gsub(".csv", replacement = "", temp)

# Catarrhini
for(i in 1:length(species)){
  genus = strsplit(species[i], "_")[[1]][1]
  sp = strsplit(species[i], "_")[[1]][2]
  
  species_subset = msrs_catarrhini[which(msrs_catarrhini$GENUS == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES == sp), ]
  
  if (nrow(species_subset) < 40) {
    print(i)
    next
  }
  
  # separate M and F
  sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]
  
  #if (nrow(sub_sexes_f) < 40 || nrow(sub_sexes_m) < 40) {
  #  print(i)
  #  next
  #}
  
  #fit = lm(log(as.matrix(species_subset[,49:87])) ~ species_subset$SEX)
  # removing PTTSP and aproxx lm
  species_subset <- subset(species_subset, select = -PTTSP)
  
  fit = lm(log(as.matrix(species_subset[,49:86])) ~ species_subset$SEX)
  cov.matrix = CalculateMatrix(fit)
  
  pasta_destino <- "~/Dropbox/Doc/Code/evowm/R/Outputs/log/" # mudar conforme necessario
  arquivo <- paste0(pasta_destino, genus, "_", sp, ".csv")
  
  write.table(cov.matrix, file = arquivo, col.names = NA)
}

# Platyrrhini
setwd("~/Dropbox/Doc/Data/p_vcv_gabriel/")
temp = list.files(pattern = "*.csv")
species  = gsub(".csv", replacement = "", temp)

# 
for(i in 1:length(species)){
  genus = strsplit(species[i], "_")[[1]][1]
  sp = strsplit(species[i], "_")[[1]][2]
  
  species_subset = msrs_platyrrhini[which(msrs_platyrrhini$GENUS. == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES. == sp), ]

  if (nrow(species_subset) < 40) {
    print(i)
    next
  }
  
  # separate M and F
  sub_sexes_m = species_subset[which(species_subset$SEX4. == "M"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX4. == "F"), ]
  
  #if (nrow(sub_sexes_f) < 20 || nrow(sub_sexes_m) < 20) {
  #  print(i)
  #  next
  #}
  species_subset <- subset(species_subset, select = -V20)
  
  fit = lm(log(as.matrix(species_subset[,23:60])) ~ species_subset$SEX4.)
  cov.matrix = CalculateMatrix(fit)
  
  pasta_destino <- "~/Dropbox/Doc/Code/evowm/R/Outputs/log/" # mudar conforme necessario
  arquivo <- paste0(pasta_destino, genus, "_", sp, ".csv")
  
  write.table(cov.matrix, file = arquivo, col.names = NA)
}
