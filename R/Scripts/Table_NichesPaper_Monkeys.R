if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(openxlsx)){install.packages("openxlsx"); library(openxlsx)}
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(Matrix)){install.packages("Matrix"); library(Matrix)}
if(!require(matrixcalc)){install.packages("matrixcalc"); library(matrixcalc)}

dados = readRDS(file = "~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

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

Respondability <- function (cov.matrix, beta.mat = NULL, iterations = 1000) {
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  respostas = apply(cov.matrix %*% beta.mat, 2, norma)
  respostas_normal = mean(respostas / tamanho_cranio)
  respostas_pure = mean(respostas)
  icv = sd(respostas) / mean(respostas)
  icv_normal = icv / tamanho_cranio
  
  lista = list(respostas, respostas_normal, icv, icv_normal, respostas_pure)
  return(lista)
}

Flexibility <- function (cov.matrix, beta.mat = NULL, iterations = 1000){
  num.traits <- dim (cov.matrix) [1]
  if(is.null(beta.mat)){
    beta.mat <- array (rnorm (num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply (beta.mat, 2, Normalize)
  }
  cor_v = vector()
  for(k in 1:iterations){
    resposta = cov.matrix %*% beta.mat[,k]
    resposta = Normalize(resposta)
    cor_v[k] = abs(corVector(resposta, beta.mat[,k]))
  }
  
  med_respo = mean(cor_v)
  icv = sd(cor_v) / mean(cor_v)
  
  lista = list(cor_v, med_respo, icv)
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
  respostas_normal = mean(respostas / tamanho_cranio)
  respostas_pure = mean(respostas)
  icv = sd(respostas) / mean(respostas)
  icv_normal = icv / tamanho_cranio
  
  lista = list(respostas, respostas_normal, icv, icv_normal, respostas_pure)
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
  respostas_normal = mean(respostas / tamanho_cranio)
  respostas_pure = mean(respostas)
  icv = sd(respostas) / mean(respostas)
  icv_normal = icv / tamanho_cranio
  
  lista = list(respostas, respostas_normal, icv, icv_normal, respostas_pure)
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
  icv = sd(respostas) / mean(respostas)
  
  lista = list(respostas, med_respo, icv)
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
  icv = sd(cor_v) / mean(cor_v)
  
  lista = list(resposta, med_respo, icv)
  return(lista)
}

##################################################################################
### Catarrhini
##################################################################################
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")
msrs[,49:87] = apply(msrs[,49:87], 2, as.numeric)

# new column with species names
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


# read all vcv matrices
setwd("~/Dropbox/Doc/Data/p_vcv_gabriel/catarrhini")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, dec = ",", header = FALSE)
names(vcv)  = gsub(".csv", replacement= "", temp)

# check matriz singularity
check = vector()
for(i in 1:length(names(vcv))){
  check[i] = is.singular.matrix(as.matrix(vcv[[i]]))
}

# check negative autovalues
check = vector()
for(i in 1:length(names(vcv))){
  autovals = eigen(as.matrix(vcv[[i]]))$values
  check[i] = length(which(autovals < 0)) >= 1
}

# check mirrored matrix equal
check = vector()
for(i in 1:length(names(vcv))){
  m = as.matrix(vcv[[i]])
  check[i] = all(m == t(m))
}

lista_subs = c("nenhum_para_ver") 

data.final = data.frame()
for(i in 1:length(temp)){
  print(i)
  # escolher sp
  genus = str_split_1(names(vcv)[[i]], "_")[1]
  sp = str_split_1(names(vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES == sp), ]
  
  # check and considering ssp
  test = sp %in% lista_subs
  if(test == TRUE){
    n_ssp = unique(species_subset$SUBSPECIES)
    
    for(k in 1:length(n_ssp)){
      ssp = str_split_1(names(vcv)[[i]], "_")[3] 
      subspecies_subset = species_subset[which(species_subset$SUBSPECIES == n_ssp[k]), ]
      
      subsp_sexes_m = subspecies_subset[which(subspecies_subset$SEX == "male"), ]
      subsp_sexes_f = subspecies_subset[which(subspecies_subset$SEX == "female"), ]
      
      subsp_med_geom = (geomean(as.matrix(subsp_sexes_m[,49:87])) + geomean(as.matrix(subsp_sexes_f[,49:87]))) / 2
      
      mean_machos = mean(geomean(as.matrix(subsp_sexes_m[,49:87])))
      mean_femeas = mean(geomean(as.matrix(subsp_sexes_f[,49:87])))
      
      # calcs evolutivos
      tamanho_cranio = sum(apply(subspecies_subset[,c(50, 55, 58, 83)] , 2, geomean))
      
      # simetrize matrix
      colnames(vcv[[i]]) = NULL
      rownames(vcv[[i]]) = NULL
      
      options(warn=-1)
      evolv = Evolvability(as.matrix(vcv[[i]]), iterations = 10000)
      evolv_n = evolv[2]
      evolv_ICV_n = evolv[3]
      evolv_pure = evolv[5]
      
      cond.evolv = ConditionalEvolvability(as.matrix(vcv[[i]]), iterations = 10000)
      cond_evolv_n = cond.evolv[2]
      cond_evolv_ICV_n = cond.evolv[3]
      cond_evolv_pure = cond.evolv[5]
      
      respond = Respondability(as.matrix(vcv[[i]]), iterations = 10000)
      respond_n = respond[2]
      respond_ICV_n = respond[3]
      respond_pure = respond[5]
      
      autonomy = Autonomy(as.matrix(vcv[[i]]), iterations = 10000)
      auto = autonomy[2]
      auto_ICV = autonomy[3]
      
      flexibility = Flexibility(as.matrix(vcv[[i]]), iterations = 10000)
      flex = flexibility[2]
      flex_ICV = flexibility[3]
      
      bias = Bias(as.matrix(vcv[[i]]), iterations = 10000)
      bias_abs = bias[2]
      bias_abs_ICV = bias[3]
      options(warn=0)
      
      # mudar para subespecies
      var_by_mean_geo = sqrt(sum(diag(as.matrix(vcv[[i]])))) / subsp_med_geom
      var_by_skull_size = sqrt(sum(diag(as.matrix(vcv[[i]])))) / tamanho_cranio
      sqrt_trace = sqrt(sum(diag(as.matrix(vcv[[i]]))))
      options(warn=-1)
      r2 = CalcR2(as.matrix(vcv[[i]]))
      options(warn=0)
      size_geomean = subsp_med_geom
      size_linham = tamanho_cranio
      traço_matriz = mean(sum(diag(as.matrix(vcv[[i]]))))
      icv = sd(eigen(vcv[[i]])$values) / mean(eigen(vcv[[i]])$values)
      cv_medio_geral = mean(apply(subspecies_subset[,49:87], 2, cv))
      cv_medio_males = mean(apply(subsp_sexes_m[,49:87], 2, cv))
      cv_medio_females = mean(apply(subsp_sexes_f[,49:87], 2, cv))
      cv_geomean_geral = cv(apply(subspecies_subset[,49:87], 1, geomean))
      cv_geomean_males = cv(apply(subsp_sexes_m[,49:87], 1, geomean))
      cv_geomean_females = cv(apply(subsp_sexes_f[,49:87], 1, geomean))
      respondability_normalized = unlist(respond_pure)
      flexibility = unlist(flex)
      evolvability_normalized = unlist(evolv_pure)
      teste_media_autovalues_evolv = mean(eigen(vcv[[i]])$values)
      cond_evolvability_normalized = unlist(cond_evolv_pure)
      autonomy = unlist(auto) 
      cv_f = unlist(flex_ICV)
      cv_e_norm = unlist(evolv_ICV_n)
      cv_c_norm = unlist(cond_evolv_ICV_n)
      cv_r_norm = unlist(respond_ICV_n)
      cv_a = unlist(auto_ICV)
      N_males = nrow(subsp_sexes_m)
      N_females = nrow(subsp_sexes_f)
      N_total = nrow(subspecies_subset)
      
      to_data = c(paste(names(vcv)[[i]], n_ssp[k], sep = "_"), var_by_mean_geo, var_by_skull_size, sqrt_trace, r2, size_geomean, size_linham, traço_matriz,
                  icv, cv_medio_geral, cv_medio_males, cv_medio_females, cv_geomean_geral, cv_geomean_males, cv_geomean_females, respondability_normalized,
                  flexibility, evolvability_normalized, teste_media_autovalues_evolv, cond_evolvability_normalized, autonomy, bias_abs, 
                  bias_abs_ICV, cv_f, cv_e_norm, cv_c_norm, cv_r_norm, cv_a, N_males, N_females, N_total, mean_machos, mean_femeas)
      
      data.final = rbind(data.final, to_data)
    }
  }
  else{
    sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
    sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]
    
    med_geom = (geomean(as.matrix(sub_sexes_m[,49:87])) + geomean(as.matrix(sub_sexes_f[,49:87]))) / 2
    
    mean_machos = mean(geomean(as.matrix(sub_sexes_m[,49:87])))
    mean_femeas = mean(geomean(as.matrix(sub_sexes_f[,49:87])))
    
    # calcs evolutivos
    tamanho_cranio = sum(apply(species_subset[,c(50, 55, 58, 83)] , 2, geomean))

    # simetrize matrix
    colnames(vcv[[i]]) = NULL
    rownames(vcv[[i]]) = NULL
    
    options(warn=-1)
    evolv = Evolvability(as.matrix(vcv[[i]]), iterations = 10000)
    evolv_n = evolv[2]
    evolv_ICV_n = evolv[3]
    evolv_pure = evolv[5]
    
    cond.evolv = ConditionalEvolvability(as.matrix(vcv[[i]]), iterations = 10000)
    cond_evolv_n = cond.evolv[2]
    cond_evolv_ICV_n = cond.evolv[3]
    cond_evolv_pure = cond.evolv[5]
    
    respond = Respondability(as.matrix(vcv[[i]]), iterations = 10000)
    respond_n = respond[2]
    respond_ICV_n = respond[3]
    respond_pure = respond[5]
    
    autonomy = Autonomy(as.matrix(vcv[[i]]), iterations = 10000)
    auto = autonomy[2]
    auto_ICV = autonomy[3]
    
    flexibility = Flexibility(as.matrix(vcv[[i]]), iterations = 10000)
    flex = flexibility[2]
    flex_ICV = flexibility[3]
    
    bias = Bias(as.matrix(vcv[[i]]), iterations = 10000)
    bias_abs = bias[2]
    bias_abs_ICV = bias[3]
    options(warn=0)
    
    var_by_mean_geo = sqrt(sum(diag(as.matrix(vcv[[i]])))) / med_geom
    var_by_skull_size = sqrt(sum(diag(as.matrix(vcv[[i]])))) / tamanho_cranio
    sqrt_trace = sqrt(sum(diag(as.matrix(vcv[[i]]))))
    options(warn=-1)
    r2 = CalcR2(as.matrix(vcv[[i]]))
    options(warn=0)
    size_geomean = med_geom
    size_linham = tamanho_cranio
    traço_matriz = mean(sum(diag(as.matrix(vcv[[i]]))))
    icv = sd(eigen(vcv[[i]])$values) / mean(eigen(vcv[[i]])$values)
    cv_medio_geral = mean(apply(species_subset[,49:87], 2, cv))
    cv_medio_males = mean(apply(sub_sexes_m[,49:87], 2, cv))
    cv_medio_females = mean(apply(sub_sexes_f[,49:87], 2, cv))
    cv_geomean_geral = cv(apply(species_subset[,49:87], 1, geomean))
    cv_geomean_males = cv(apply(sub_sexes_m[,49:87], 1, geomean))
    cv_geomean_females = cv(apply(sub_sexes_f[,49:87], 1, geomean))
    respondability_normalized = unlist(respond_pure)
    flexibility = unlist(flex)
    evolvability_normalized = unlist(evolv_pure)
    teste_media_autovalues_evolv = mean(eigen(vcv[[i]])$values)
    cond_evolvability_normalized = unlist(cond_evolv_pure)
    autonomy = unlist(auto) 
    cv_f = unlist(flex_ICV)
    cv_e_norm = unlist(evolv_ICV_n)
    cv_c_norm = unlist(cond_evolv_ICV_n)
    cv_r_norm = unlist(respond_ICV_n)
    cv_a = unlist(auto_ICV)
    N_males = nrow(sub_sexes_m)
    N_females = nrow(sub_sexes_f)
    N_total = nrow(species_subset)
    
    to_data = c(names(vcv)[[i]], var_by_mean_geo, var_by_skull_size, sqrt_trace, r2, size_geomean, size_linham, traço_matriz,
                icv, cv_medio_geral, cv_medio_males, cv_medio_females, cv_geomean_geral, cv_geomean_males, cv_geomean_females, respondability_normalized,
                flexibility, evolvability_normalized, teste_media_autovalues_evolv, cond_evolvability_normalized, autonomy, bias_abs, 
                bias_abs_ICV, cv_f, cv_e_norm, cv_c_norm, cv_r_norm, cv_a, N_males, N_females, N_total, mean_machos, mean_femeas)
    
    data.final = rbind(data.final, to_data)
  }
}

colnames(data.final) = c("Especie", "var_by_mean_geo", "var_by_skull_size", "sqrt_trace", "r2", "size_geomean", "size_linham", "traço_matriz",
                         "icv",  "cv_medio_geral", "cv_medio_males", "cv_medio_females", "cv_geomean_geral", "cv_geomean_males", "cv_geomean_females", "respondability_normalized",
                         "flexibility", "evolvability_normalized", "teste_media_autovalues_evolv", "cond_evolvability_normalized", 
                         "autonomy", "bias_abs", "bias_abs_CV", "cv_f", "cv_e_normalized", "cv_c_normalized", "cv_r_normalized", 
                         "cv_a", "N_males", "N_females", "N_total", "media_geom_machos", "media_geom_femeas")


write.xlsx(data.final, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/Table_NichesPaper_Catarrhini.xlsx")

##################################################################################
### Platyrrhini
##################################################################################
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini.csv", dec = ",", sep = ",")
msrs[,23:61] = apply(msrs[,23:61], 2, as.numeric)

# new column with species names
species = vector()
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUB.[i])){
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], msrs$SUB.[i], sep="_")
  }
}

species = unique(species)

index = which(apply(msrs[, 23:61], 1, function(x) any(is.na(x) | (x == "" & !is.numeric(x)))))
msrs = msrs[-index,]

# remove uncertain sex
msrs$SEX4.[msrs$SEX4. == ""] = NA
msrs$SEX4.[msrs$SEX4. == " "] = NA

msrs = msrs[complete.cases(msrs$SEX4.), ]

# read all vcv matrices
setwd("~/Dropbox/Doc/Data/p_vcv_gabriel")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, dec = ",", header = FALSE)
names(vcv)  = gsub(".csv", replacement= "", temp)

# check matriz singularity
check = vector()
for(i in 1:length(names(vcv))){
  check[i] = is.singular.matrix(as.matrix(vcv[[i]]))
}

# check negative autovalues
check = vector()
for(i in 1:length(names(vcv))){
  autovals = eigen(as.matrix(vcv[[i]]))$values
  check[i] = length(which(autovals < 0)) >= 1
}

# check mirrored matrix equal
check = vector()
for(i in 1:length(names(vcv))){
  m = as.matrix(vcv[[i]])
  check[i] = all(m == t(m))
}

lista_subs = c("belzebul", "palliata", "personatus", "torquatus", "libidinosus", "nigritus", 
               "irrorata", "pithecia", "midas", "mystax", "nigricollis")

data.final = data.frame()
for(i in 1:length(temp)){
  print(i)
  # escolher sp
  genus = str_split_1(names(vcv)[[i]], "_")[1]
  sp = str_split_1(names(vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS. == genus), ] 
  species_subset = species_subset[which(species_subset$SPECIES. == sp), ] 
  
  # check and considering ssp
  test = sp %in% lista_subs
  if(test == TRUE){
    n_ssp = unique(species_subset$SUB.)
    
    for(k in 1:length(n_ssp)){
      ssp = str_split_1(names(vcv)[[i]], "_")[3] 
      subspecies_subset = species_subset[which(species_subset$SUB. == n_ssp[k]), ]
      
      subsp_sexes_m = subspecies_subset[which(subspecies_subset$SEX4. == "M"), ]
      subsp_sexes_f = subspecies_subset[which(subspecies_subset$SEX4. == "F"), ]
      
      subsp_med_geom = (geomean(as.matrix(subsp_sexes_m[,23:61])) + geomean(as.matrix(subsp_sexes_f[,23:61]))) / 2
      
      media_machos = mean(geomean(as.matrix(subsp_sexes_m[,23:61])))
      media_femeas = mean(geomean(as.matrix(subsp_sexes_f[,23:61])))
      
      # calcs evolutivos
      tamanho_cranio = sum(apply(subspecies_subset[,c(24, 29, 32, 57)] , 2, geomean))
      
      # simetrize matrix
      colnames(vcv[[i]]) = NULL
      rownames(vcv[[i]]) = NULL
      
      options(warn=-1)
      evolv = Evolvability(as.matrix(vcv[[i]]), iterations = 10000)
      evolv_n = evolv[2]
      evolv_ICV_n = evolv[3]
      evolv_pure = evolv[5]
      
      cond.evolv = ConditionalEvolvability(as.matrix(vcv[[i]]), iterations = 10000)
      cond_evolv_n = cond.evolv[2]
      cond_evolv_ICV_n = cond.evolv[3]
      cond_evolv_pure = cond.evolv[5]
      
      respond = Respondability(as.matrix(vcv[[i]]), iterations = 10000)
      respond_n = respond[2]
      respond_ICV_n = respond[3]
      respond_pure = respond[5]
      
      autonomy = Autonomy(as.matrix(vcv[[i]]), iterations = 10000)
      auto = autonomy[2]
      auto_ICV = autonomy[3]
      
      flexibility = Flexibility(as.matrix(vcv[[i]]), iterations = 10000)
      flex = flexibility[2]
      flex_ICV = flexibility[3]
      
      bias = Bias(as.matrix(vcv[[i]]), iterations = 10000)
      bias_abs = bias[2]
      bias_abs_ICV = bias[3]
      options(warn=0)
      
      # mudar para subespecies
      var_by_mean_geo = sqrt(sum(diag(as.matrix(vcv[[i]])))) / subsp_med_geom
      var_by_skull_size = sqrt(sum(diag(as.matrix(vcv[[i]])))) / tamanho_cranio
      sqrt_trace = sqrt(sum(diag(as.matrix(vcv[[i]]))))
      options(warn=-1)
      r2 = CalcR2(as.matrix(vcv[[i]]))
      options(warn=0)
      size_geomean = subsp_med_geom
      size_linham = tamanho_cranio
      traço_matriz = mean(sum(diag(as.matrix(vcv[[i]]))))
      icv = sd(eigen(vcv[[i]])$values) / mean(eigen(vcv[[i]])$values)
      cv_medio_geral = mean(apply(subspecies_subset[,23:61] , 2, cv))
      cv_medio_males = mean(apply(subsp_sexes_m[,23:61] , 2, cv))
      cv_medio_females = mean(apply(subsp_sexes_f[,23:61] , 2, cv))
      cv_geomean_geral = cv(apply(subspecies_subset[,23:61] , 1, geomean))
      cv_geomean_males = cv(apply(subsp_sexes_m[,23:61] , 1, geomean))
      cv_geomean_females = cv(apply(subsp_sexes_f[,23:61] , 1, geomean))
      respondability_normalized = unlist(respond_pure)
      flexibility = unlist(flex)
      evolvability_normalized = unlist(evolv_pure)
      teste_media_autovalues_evolv = mean(eigen(vcv[[i]])$values)
      cond_evolvability_normalized = unlist(cond_evolv_pure)
      autonomy = unlist(auto) 
      cv_f = unlist(flex_ICV)
      cv_e_norm = unlist(evolv_ICV_n)
      cv_c_norm = unlist(cond_evolv_ICV_n)
      cv_r_norm = unlist(respond_ICV_n)
      cv_a = unlist(auto_ICV)
      N_males = nrow(subsp_sexes_m)
      N_females = nrow(subsp_sexes_f)
      N_total = nrow(subspecies_subset)
      
      to_data = c(paste(names(vcv)[[i]], n_ssp[k], sep = "_"), var_by_mean_geo, var_by_skull_size, sqrt_trace, r2, size_geomean, size_linham, traço_matriz,
                  icv, cv_medio_geral, cv_medio_males, cv_medio_females, cv_geomean_geral, cv_geomean_males, cv_geomean_females, respondability_normalized,
                  flexibility, evolvability_normalized, teste_media_autovalues_evolv, cond_evolvability_normalized, autonomy, bias_abs, 
                  bias_abs_ICV, cv_f, cv_e_norm, cv_c_norm, cv_r_norm, cv_a, N_males, N_females, N_total, media_machos, media_femeas)
      
      data.final = rbind(data.final, to_data)
    }
  }
  else{
    sub_sexes_m = species_subset[which(species_subset$SEX4. == "M"), ]
    sub_sexes_f = species_subset[which(species_subset$SEX4. == "F"), ]
  
    med_geom = (geomean(as.matrix(sub_sexes_m[,23:61])) + geomean(as.matrix(sub_sexes_f[,23:61]))) / 2
    
    media_machos = mean(geomean(as.matrix(sub_sexes_m[,23:61])))
    media_femeas = mean(geomean(as.matrix(sub_sexes_f[,23:61])))
    
    # calcs evolutivos
    tamanho_cranio = sum(apply(species_subset[,c(24, 29, 32, 57)] , 2, geomean)) 
  
    # simetrize matrix
    colnames(vcv[[i]]) = NULL
    rownames(vcv[[i]]) = NULL
  
    options(warn=-1)
    evolv = Evolvability(as.matrix(vcv[[i]]), iterations = 10000)
    evolv_n = evolv[2]
    evolv_ICV_n = evolv[3]
    evolv_pure = evolv[5]
  
    cond.evolv = ConditionalEvolvability(as.matrix(vcv[[i]]), iterations = 10000)
    cond_evolv_n = cond.evolv[2]
    cond_evolv_ICV_n = cond.evolv[3]
    cond_evolv_pure = cond.evolv[5]
  
    respond = Respondability(as.matrix(vcv[[i]]), iterations = 10000)
    respond_n = respond[2]
    respond_ICV_n = respond[3]
    respond_pure = respond[5]
  
    autonomy = Autonomy(as.matrix(vcv[[i]]), iterations = 10000)
    auto = autonomy[2]
    auto_ICV = autonomy[3]
  
    flexibility = Flexibility(as.matrix(vcv[[i]]), iterations = 10000)
    flex = flexibility[2]
    flex_ICV = flexibility[3]
  
    bias = Bias(as.matrix(vcv[[i]]), iterations = 10000)
    bias_abs = bias[2]
    bias_abs_ICV = bias[3]
    options(warn=0)
  
    var_by_mean_geo = sqrt(sum(diag(as.matrix(vcv[[i]])))) / med_geom
    var_by_skull_size = sqrt(sum(diag(as.matrix(vcv[[i]])))) / tamanho_cranio
    sqrt_trace = sqrt(sum(diag(as.matrix(vcv[[i]]))))
    options(warn=-1)
    r2 = CalcR2(as.matrix(vcv[[i]]))
    options(warn=0)
    size_geomean = med_geom
    size_linham = tamanho_cranio
    traço_matriz = mean(sum(diag(as.matrix(vcv[[i]]))))
    icv = sd(eigen(vcv[[i]])$values) / mean(eigen(vcv[[i]])$values)
    cv_medio_geral = mean(apply(species_subset[,23:61] , 2, cv))
    cv_medio_males = mean(apply(sub_sexes_m[,23:61] , 2, cv))
    cv_medio_females = mean(apply(sub_sexes_f[,23:61] , 2, cv))
    cv_geomean_geral = cv(apply(species_subset[,23:61] , 1, geomean))
    cv_geomean_males = cv(apply(sub_sexes_m[,23:61] , 1, geomean))
    cv_geomean_females = cv(apply(sub_sexes_f[,23:61] , 1, geomean))
    respondability_normalized = unlist(respond_pure)
    flexibility = unlist(flex)
    evolvability_normalized = unlist(evolv_pure)
    teste_media_autovalues_evolv = mean(eigen(vcv[[i]])$values)
    cond_evolvability_normalized = unlist(cond_evolv_pure)
    autonomy = unlist(auto) 
    cv_f = unlist(flex_ICV)
    cv_e_norm = unlist(evolv_ICV_n)
    cv_c_norm = unlist(cond_evolv_ICV_n)
    cv_r_norm = unlist(respond_ICV_n)
    cv_a = unlist(auto_ICV)
    N_males = nrow(sub_sexes_m)
    N_females = nrow(sub_sexes_f)
    N_total = nrow(species_subset)
    
    to_data = c(names(vcv)[[i]], var_by_mean_geo, var_by_skull_size, sqrt_trace, r2, size_geomean, size_linham, traço_matriz,
                icv, cv_medio_geral, cv_medio_males, cv_medio_females, cv_geomean_geral, cv_geomean_males, cv_geomean_females, respondability_normalized,
                flexibility, evolvability_normalized, teste_media_autovalues_evolv, cond_evolvability_normalized, autonomy, bias_abs, 
                bias_abs_ICV, cv_f, cv_e_norm, cv_c_norm, cv_r_norm, cv_a, N_males, N_females, N_total, media_machos, media_femeas)
  
    data.final = rbind(data.final, to_data)
  }
}

colnames(data.final) = c("Especie", "var_by_mean_geo", "var_by_skull_size", "sqrt_trace", "r2", "size_geomean", "size_linham", "traço_matriz",
                         "icv",  "cv_medio_geral", "cv_medio_males", "cv_medio_females", "cv_geomean_geral", "cv_geomean_males", "cv_geomean_females", "respondability_normalized",
                         "flexibility", "evolvability_normalized", "teste_media_autovalues_evolv", "cond_evolvability_normalized", 
                         "autonomy", "bias_abs", "bias_abs_CV", "cv_f", "cv_e_normalized", "cv_c_normalized", "cv_r_normalized", 
                         "cv_a", "N_males", "N_females", "N_total", "media_geom_machos", "media_geom_femeas")


write.xlsx(data.final, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/Table_NichesPaper_Platyrrhini.xlsx")
