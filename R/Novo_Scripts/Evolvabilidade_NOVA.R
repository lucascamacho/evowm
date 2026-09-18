##################################################################################
### Catarrhini
##################################################################################
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(openxlsx)){install.packages("openxlsx"); library(openxlsx)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(Matrix)){install.packages("Matrix"); library(Matrix)}
if(!require(matrixcalc)){install.packages("matrixcalc"); library(matrixcalc)}

#msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")
msrs = read.xlsx("~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini_denovo.xlsx", colNames = TRUE, check.names = TRUE)
msrs[,49:87] = apply(msrs[,49:87], 2, as.numeric)

# species names
species = vector()
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUBSPECIES[i])){
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], msrs$SUBSPECIES[i], sep="_")
  }
}

species = unique(species)

# check names and remove doubts
msrs$SEX[which(msrs$SEX == "?female")] = "female"
msrs$SEX[which(msrs$SEX == "?male")] = "male"

# remove uncertain sex
msrs$SEX[msrs$SEX == "0"] = NA
msrs$SEX[msrs$SEX == ""] = NA
msrs$SEX[msrs$SEX == "sexo"] = NA

msrs = msrs[complete.cases(msrs$SEX), ]

# evo functions
Normalize <- function(x){return(x/Norm(x))}
Norm <- function(x){return(sqrt(sum(x*x)))}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# Definindo parâmetros
#n_linhas <- 39
#n_colunas <- 10000
#media <- 0
#desvio_padrao <- 1

# Gerando a matriz com rnorm
#matriz <- matrix(rnorm(n_linhas * n_colunas, mean = media, sd = desvio_padrao), 
#                 nrow = n_linhas, ncol = n_colunas)
#beta = apply (matriz, 2, Normalize)

Evolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  #respostas = sum((beta.mat[,1] %*% cov.matrix) * beta.mat[,1])
  #respostas = diag(t(beta.mat) %*% cov.matrix %*% beta.mat)
  respostas = diag(t(beta.mat) %*% cov.matrix %*% beta.mat) / med_geom
  respostas_normal = mean(respostas / tamanho_cranio)
  respostas_pure = mean(respostas)
  icv = sd(respostas) / mean(respostas)
  icv_normal = icv / tamanho_cranio
  
  lista = list(respostas, respostas_normal, icv, icv_normal, respostas_pure)
  return(lista)
}

# read all vcv matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Novo_Output/p_vcv_gabriel/catarrhini")
temp = list.files(pattern = "*.csv")
nolog_vcv = lapply(temp, read.csv, dec = ",", header = FALSE)
names(nolog_vcv)  = gsub(".csv", replacement= "", temp)

data.final = data.frame()
for(i in 1:length(temp)){
  print(i)
  # escolher sp
  genus = str_split_1(names(nolog_vcv)[[i]], "_")[1]
  sp = str_split_1(names(nolog_vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES == sp), ]
  
  sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]
  
  # calcs evolutivos
  tamanho_cranio = sum(apply(species_subset[,c(50, 55, 58, 83)] , 2, geomean))
  med_geom = (geomean(as.matrix(sub_sexes_m[,49:87])) + geomean(as.matrix(sub_sexes_f[,49:87]))) / 2
  
  # simetrize matrix
  colnames(nolog_vcv[[i]]) = NULL
  rownames(nolog_vcv[[i]]) = NULL
  
  evolv = Evolvability(as.matrix(nolog_vcv[[i]]), iterations = 10000)
  
  hist(evolv[[1]])
  evolv_mean = mean(evolv[[1]])
  qs_5 = quantile(evolv[[1]], probs = c(0.05, 0.95))[1]
  qs_95 = quantile(evolv[[1]], probs = c(0.05, 0.95))[2]
  
  
  # dimor evolvability
  dimor = apply(as.matrix(sub_sexes_m[,49:87]), 2, geomean) - apply(as.matrix(sub_sexes_f[,49:87]), 2, geomean)
  
  # evolvability
  dimor_norm = dimor / sqrt(sum(dimor ^ 2))

  #resposta = as.matrix(nolog_vcv[[i]]) %*% as.vector(dimor_norm)
  #evolv = sum(resposta * dimor_norm)
  
  cov.matrix = as.matrix(nolog_vcv[[i]]) 
  evolvi = sum((dimor_norm %*% cov.matrix) * dimor_norm)
  
  # std evolvability
  evolvi_norm = sqrt(evolvi) / med_geom

  hist(evolv[[1]], main = paste(genus, sp))
  abline(v = qs_5, col = "red", lwd = 3)
  abline(v = qs_95, col = "red", lwd = 3)
  abline(v = evolvi_norm, col = "blue", lwd = 3)
  
  # check quantils for plot
  if (evolvi_norm > qs_95 | evolvi_norm < qs_5) {
    check = 1
  } else {
    check = 0
  }
  
  to_data = c(names(nolog_vcv)[[i]], mean(dimor), evolv_mean, qs_5, qs_95, evolvi, evolvi_norm, check)
  data.final = rbind(data.final, to_data)
}

colnames(data.final) = c("species", "media_dimorfismo", "evolv_nula_media", "qs_5", "qs_95", "evolv_dimor", "evolv_dimor_norm", "check")

##################################################################################
### Ṕlatyrrhini
##################################################################################
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(openxlsx)){install.packages("openxlsx"); library(openxlsx)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(Matrix)){install.packages("Matrix"); library(Matrix)}
if(!require(matrixcalc)){install.packages("matrixcalc"); library(matrixcalc)}

#msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini_2.csv", dec = ",", sep = ",")
msrs = read.xlsx("~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini_denovo.xlsx", colNames = TRUE, check.names = TRUE)
msrs[,23:61] = apply(msrs[,23:61], 2, as.numeric)
#msrs$V1 = as.numeric(msrs$V1)

# new column with species names
species = vector()
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUB.[i])){
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], msrs$SUB.[i], sep="_")
  }
}

#msrs = cbind(species, msrs)
species = unique(species)

index = which(apply(msrs[,23:61], 1, anyNA) == TRUE)
msrs = msrs[-index,]

# remove uncertain sex
msrs$SEX4.[msrs$SEX4. == ""] = NA
msrs$SEX4.[msrs$SEX4. == " "] = NA

msrs = msrs[complete.cases(msrs$SEX4.), ]

# evo functions
Normalize <- function(x){return(x/Norm(x))}
Norm <- function(x){return(sqrt(sum(x*x)))}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# Definindo parâmetros
n_linhas <- 39
n_colunas <- 10000
media <- 0
desvio_padrao <- 1

# Gerando a matriz com rnorm
matriz <- matrix(rnorm(n_linhas * n_colunas, mean = media, sd = desvio_padrao), 
                 nrow = n_linhas, ncol = n_colunas)
beta = apply (matriz, 2, Normalize)

Evolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  #respostas = sum((beta.mat[,1] %*% cov.matrix) * beta.mat[,1])
  respostas = diag(t(beta.mat) %*% cov.matrix %*% beta.mat)
  respostas_normal = mean(respostas / tamanho_cranio)
  respostas_pure = mean(respostas)
  icv = sd(respostas) / mean(respostas)
  icv_normal = icv / tamanho_cranio
  
  lista = list(respostas, respostas_normal, icv, icv_normal, respostas_pure)
  return(lista)
}

# read all vcv matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Novo_Output/p_vcv_gabriel")
temp = list.files(pattern = "*.csv")
nolog_vcv = lapply(temp, read.csv, dec = ",", header = FALSE)
names(nolog_vcv)  = gsub(".csv", replacement= "", temp)

data.final = data.frame()
for(i in 1:length(temp)){
  print(i)
  # escolher sp
  genus = str_split_1(names(nolog_vcv)[[i]], "_")[1]
  sp = str_split_1(names(nolog_vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS. == genus), ] 
  species_subset = species_subset[which(species_subset$SPECIES. == sp), ] 
  
  sub_sexes_m = species_subset[which(species_subset$SEX4. == "M"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX4. == "F"), ]
  
  # calcs evolutivos
  tamanho_cranio = sum(apply(species_subset[,c(24, 29, 32, 57)] , 2, geomean)) 
  med_geom = (geomean(as.matrix(sub_sexes_m[,23:61])) + geomean(as.matrix(sub_sexes_f[,23:61]))) / 2
  
  # simetrize matrix
  colnames(nolog_vcv[[i]]) = NULL
  rownames(nolog_vcv[[i]]) = NULL
  
  evolv = Evolvability(as.matrix(nolog_vcv[[i]]), beta.mat = beta, iterations = 10000)
  hist(evolv[[1]])
  evolv_mean = mean(evolv[[1]])
  qs_5 = quantile(evolv[[1]], probs = c(0.05, 0.95))[1]
  qs_95 = quantile(evolv[[1]], probs = c(0.05, 0.95))[2]
  
  # dimor evolvability
  dimor = apply(as.matrix(sub_sexes_m[,23:61]), 2, geomean) - apply(as.matrix(sub_sexes_f[,23:61]), 2, geomean)
  
  # evolvability
  dimor_norm = dimor / sqrt(sum(dimor ^ 2))
  
  #resposta = as.matrix(nolog_vcv[[i]]) %*% as.vector(dimor_norm)
  
  #evolv = sum(resposta * dimor_norm)
  cov.matrix = as.matrix(nolog_vcv[[i]]) 
  evolv = sum((dimor_norm %*% cov.matrix) * dimor_norm)
  
  # std evolvability
  evolv_norm = sqrt(evolv) / med_geom
  
  # check quantils for plot
  if (evolv > qs_95 | evolv < qs_5) {
    check = 1
  } else {
    check = 0
  }
  
  to_data = c(names(nolog_vcv)[[i]], mean(dimor), evolv_mean, qs_5, qs_95, evolv, evolv_norm, check)
  data.final = rbind(data.final, to_data)
}

colnames(data.final) = c("species", "media_dimorfismo", "evolv_nula_media", "qs_5", "qs_95", "evolv_dimor", "evolv_dimor_norm", "check")

