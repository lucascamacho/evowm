# set WD and functions
setwd("~/Dropbox/Doc/Output/")

# load packages
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

pars = read.table("parsimonia_quadrada.csv", sep = ",", header = TRUE)
ml = read.table("ml.txt", sep = ",", header = TRUE)

corela = vector()
for(i in 2:40){
  corela[i] = corVector(pars[,i], ml[,i])
}

corela
mean(corela[2:40])
