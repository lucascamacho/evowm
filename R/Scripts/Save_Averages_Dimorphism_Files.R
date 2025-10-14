# # set WD and functions
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

# load packages
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# load data to correlation
medidas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

# dimorphism
medias_machos = data.frame()
medias_femeas = data.frame()
for(i in 1:length(medidas$Species)){
  medias_machos = rbind(medias_machos, medidas$ByTrait_Averages[[i]]$Machos)
  medias_femeas = rbind(medias_femeas, medidas$ByTrait_Averages[[i]]$Fêmeas)

}

# insert species names and traits names
medias_machos = cbind(medidas$Species, medias_machos)
medias_femeas = cbind(medidas$Species, medias_femeas)

colnames(medias_machos)[2:40] = names(medidas$ByTrait_Averages$Cercocebus_torquatus$Machos)
colnames(medias_femeas)[2:40] = names(medidas$ByTrait_Averages$Cercocebus_torquatus$Fêmeas)

colnames(medias_machos)[1] = "Species" 
colnames(medias_femeas)[1] = "Species" 

# save averages tables
#write.csv(medias_machos, file = "avg_males.csv")
#write.csv(medias_femeas, file = "avg_females.csv")

# save dimorphism table
dimor = as.matrix(medias_machos[,2:40]) - as.matrix(medias_femeas[,2:40])
dimor = cbind(medidas$Species, dimor)
colnames(dimor)[1] = "Species"

#write.csv(dimor, file = "avg_dimorphism.csv")