##################################################################################
### Catarrhini
##################################################################################


if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(openxlsx)){install.packages("openxlsx"); library(openxlsx)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(Matrix)){install.packages("Matrix"); library(Matrix)}

dados = readRDS(file = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Averages_PCS_Extant_Species.RData")

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

cv = function(x) sd(x)/mean(x)

#

Respondability <- function (cov.matrix, beta.mat = NULL, iterations = 1000) {
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations)) #sample de 39 traits de uma norm
    beta.mat <- apply (beta.mat, 2, Normalize) # normaliza todos os vetores normais
  }
  respostas = apply(cov.matrix %*% beta.mat, 2, norma) # multiplica pela matriz VCV e pega a norma
  med_respo = mean(respostas)
  icv_respo = sd(respostas) / mean(respostas)
  
  lista = list(respostas, med_respo, icv_respo)
  return(lista)
}

Flexibility <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  #Cb <- apply(cov.matrix %*% beta.mat, 2, Normalize)
  #respostas = diag(t (beta.mat) %*% Cb)
  cor_v = vector()
  for(k in 1:iterations){
    resposta = cov.matrix %*% beta.mat[,k]
    resposta = Normalize(resposta)
    cor_v[k] = abs(corVector(resposta, beta.mat[,k]))
  }
  
  med_respo = mean(cor_v)
  icv_respo = sd(cor_v)/mean(cor_v)
  
  lista = list(cor_v, med_respo, icv_respo)
  return(lista)
}

Evolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  #respostas = sum((beta.mat[,1] %*% cov.matrix) * beta.mat[,1])
  respostas = diag(t(beta.mat) %*% cov.matrix %*% beta.mat)
  
  med_respo = mean(respostas)
  icv_respo = sd(respostas)/mean(respostas)
  
  lista = list(respostas, med_respo, icv_respo)
  return(lista)
}

ConditionalEvolvability <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  cov.matrix <- Matrix(cov.matrix)
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  cov.matrix = tryCatch({chol(cov.matrix); cov.matrix}, error = function(cond){
    warning("matrix is singular, can't compute conditional evolvability directly. Using nearPD, results could be wrong")
    cov.matrix <- nearPD(cov.matrix)[[1]]
    chol(cov.matrix)
  })
  #respostas = sum((beta.mat[,1] %*% solve(cov.matrix)) * beta.mat[,1])
  respostas =  (1/diag(t (beta.mat) %*% solve (cov.matrix, beta.mat)))
  
  med_respo = mean(respostas)
  icv_respo = sd(respostas)/mean(respostas)
  
  lista = list(respostas, med_respo, icv_respo)
  return(lista)
}

Autonomy <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  cov.matrix <- Matrix(cov.matrix)
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  cov.matrix = tryCatch({cv = chol(cov.matrix); cov.matrix}, error = function(cond){
    warning("matrix is singular, can't compute autonomy directly. Using nearPD, results could be wrong")
    cov.matrix <- nearPD(cov.matrix)[[1]]
    chol(cov.matrix)
  })
  respostas = (1/diag(t (beta.mat) %*% solve (cov.matrix, beta.mat))) / diag(t(beta.mat) %*% cov.matrix %*% beta.mat)
  
  med_respo = mean(respostas)
  icv_respo = sd(respostas)/mean(respostas)
  
  lista = list(respostas, med_respo, icv_respo)
  return(lista)
}

Bias <- function (cov.matrix, beta.mat = NULL, iterations = 1000) {
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations)) #sample de 39 traits de uma norm
    beta.mat <- apply (beta.mat, 2, Normalize) # normaliza todos os vetores normais
  }
  
  cor_v = vector()
  for(k in 1:iterations){
    resposta = cov.matrix %*% beta.mat[,k]
    resposta = Normalize(resposta)
    cor_v[k] = abs(corVector(resposta, eigen(cov.matrix)$vectors[,1]))
  }
  
  med_respo = mean(cor_v)
  icv_respo = sd(cor_v)/mean(cor_v)
  
  lista = list(resposta, med_respo, icv_respo)
  return(lista)
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

# read all nolog vcv matrices
setwd("~/Dropbox/Doc/Output/nolog_vcv")
temp = list.files(pattern = "*.txt")
# remove species which are not in the phylogeny
to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")
index = grep(paste(to_remove, collapse = "|"), temp)
temp = temp[-index]
nolog_vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(nolog_vcv)  = gsub(".txt", replacement= "", temp)

data.final = data.frame()
for(i in 1:length(dados$Species)){
  # SOMENTE NO LOG
  # geo mean
  nolog_med_geom =  geomean((dados$NoLog_ByTrait_Averages[[i]]$Machos + 
                       dados$NoLog_ByTrait_Averages[[i]]$Fêmeas) / 2)
  
  # calcs evolutivos
  options(warn=-1)
  evolv = Evolvability(as.matrix(nolog_vcv[[i]]), iterations = 10000)
  nolog_evolv = evolv[2]
  nolog_evolv_ICV = evolv[3]
  cond.evolv = ConditionalEvolvability(as.matrix(nolog_vcv[[i]]), iterations = 10000)
  nolog_cond_evolv = cond.evolv[2]
  nolog_cond_evolv_ICV = cond.evolv[3]
  respond = Respondability(as.matrix(nolog_vcv[[i]]), iterations = 10000)
  nolog_respond = respond[2]
  nolog_respond_ICV = respond[3]
  auto = Autonomy(as.matrix(nolog_vcv[[i]]), iterations = 10000)
  nolog_auto = auto[2]
  nolog_auto_ICV = auto[3]
  flex = Flexibility(as.matrix(nolog_vcv[[i]]), iterations = 10000)
  nolog_flex = flex[2]
  nolog_flex_ICV = flex[3]
  bias = Bias(as.matrix(nolog_vcv[[i]]), iterations = 10000)
  bias_abs = bias[2]
  nolog_bias_ICV = bias[3]
  options(warn=0)
  
  # escolher
  genus = str_split_1(names(nolog_vcv)[[i]], "_")[1]
  sp = str_split_1(names(nolog_vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES == sp), ]
  
  sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]
  
  # calcs evolutivos
  tamanho_cranio = sum(apply(species_subset[,c(50, 55, 58, 83)] , 2, geomean))
  var_by_mean_geo = sqrt(sum(diag(as.matrix(nolog_vcv[[i]])))) / nolog_med_geom
  var_by_skull_size = sqrt(sum(diag(as.matrix(nolog_vcv[[i]])))) / tamanho_cranio
  sqrt_trace = sqrt(sum(diag(as.matrix(nolog_vcv[[i]]))))
  options(warn=-1)
  r2 = CalcR2(as.matrix(nolog_vcv[[i]]))
  options(warn=0)
  size_geomean = nolog_med_geom
  size_linham = mean(sum(diag(as.matrix(nolog_vcv[[i]]))))
  icv = sd(eigen(nolog_vcv[[i]])$values) / mean(eigen(nolog_vcv[[i]])$values)
  cv_medio = mean(apply(species_subset[,49:87] , 2, cv))
  cv_geomean_geral = cv(apply(species_subset[,49:87] , 2, geomean))
  cv_geomean_males = cv(apply(sub_sexes_m[,49:87] , 2, geomean))
  cv_geomean_females = cv(apply(sub_sexes_f[,49:87] , 2, geomean))
  evolvability = unlist(nolog_evolv) / tamanho_cranio
  cond_evolvability = unlist(nolog_cond_evolv)
  flexibility = unlist(nolog_flex)
  respondability = unlist(nolog_respond) / tamanho_cranio
  autonomy = unlist(nolog_auto)
  icv_e_norm = unlist(nolog_evolv_ICV) / tamanho_cranio
  icv_c_norm = unlist(nolog_cond_evolv_ICV) / tamanho_cranio
  icv_r_norm = unlist(nolog_respond_ICV) / tamanho_cranio
  
    
  to_data = c(dados$Species[[i]], var_by_mean_geo, var_by_skull_size, sqrt_trace, r2, size_geomean, size_linham,
              icv, cv_medio, cv_geomean_geral, cv_geomean_males, cv_geomean_females, respondability, flexibility,
              nolog_evolv, nolog_cond_evolv, bias_abs, evolvability, cond_evolvability, autonomy, nolog_respond_ICV,
              nolog_flex_ICV, nolog_evolv_ICV, nolog_cond_evolv_ICV, icv_e_norm, icv_c_norm, icv_r_norm, nolog_auto_ICV)
    
  data.final = rbind(data.final, to_data)
}

colnames(data.final) = c("Especie", "var_by_mean_geo", "var_by_skull_size", "sqrt_trace", "r2", "size_geomean", "size_linham",
                         "icv", "cv_medio", "cv_geomean_geral", "cv_geomean_males", "cv_geomean_females", "respondability", "flex",
                         "evolvabilidade", "conditional_evolvability", "Bias_abs", "E_normalized", "C_normalized", "Autonomy", 
                         "CV_R", "CV_Flex", "CV_E", "CV_C", "CV_Enorm", "CV_C_norm", "CV_R_norm", "CV_A")

  
write.xlsx(data.final, file = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Table_paranicho_Catarrhini.xlsx")

rm(list = ls())


