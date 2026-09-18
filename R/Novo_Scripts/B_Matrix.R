# set WD and functions
setwd("~/Dropbox/Doc/")

# load packages
if(!require(stringr)){install.packages("stringr"); library(stringr)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(matrixcalc)){install.packages("matrixcalc"); library(matrixcalc)}

# Função para simetrizar a matriz
simetrizar_matriz <- function(mat) {
  return((mat + t(mat)) / 2)
}

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

# Read and unique species
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")

setwd("~/Dropbox/Doc/Code/evowm/R/Novo_Output/p_vcv_gabriel/catarrhini")
temp = list.files(pattern = "*.csv")
nolog_vcv = lapply(temp, read.csv, dec = ",", header = FALSE)
names(nolog_vcv)  = gsub(".csv", replacement= "", temp)

# new column with species names
species = vector()
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUBSPECIES[i])){
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], msrs$SUBSPECIES[i], sep="_")
  }
}

# Converter para matriz de correlação
cov_para_cor <- function(B) {
  desvios <- sqrt(diag(B))
  desvios[desvios < .Machine$double.eps] <- NA
  R <- B / outer(desvios, desvios)
  R[is.na(R)] <- 0
  return(R)
}

# Calcular matriz B corretamente
calcular_B_manova <- function(species_subset) {
  # Ajustar o modelo MANOVA
  modelo <- manova(as.matrix(species_subset[, 50:88]) ~ species_subset$SEX)
  
  # Extrair matriz B (matriz H: soma de quadrados entre grupos)
  b <- summary(modelo, test = "Wilks")$SS$'species_subset$SEX'
  
  # Garantir simetria numérica
  b <- (b + t(b)) / 2
  
  # Normalizar por (n - 1)
  return(b / (nrow(species_subset) - 1))
}


msrs = cbind(species, msrs)
species = unique(species)


###############################################################################################################


bes = list()
bes_medias = list()
for(i in 1:length(temp)){
  i = 34
  # escolher sp
  genus = str_split_1(names(nolog_vcv)[[i]], "_")[1]
  sp = str_split_1(names(nolog_vcv)[[i]], "_")[2]
  
  species_subset = msrs[which(msrs$GENUS == genus), ]
  species_subset = species_subset[which(species_subset$SPECIES == sp), ]
  
  print(nrow(species_subset))
  
  sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
  sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]
  
  bes[[i]] = x
  b_manova = calcular_B_manova(species_subset)
  
  
  corrigir_com_pca <- function(B) {
    # Decomposição em autovalores e autovetores
    eig <- eigen(B)
    
    # Zerando autovalores negativos
    eig$values[eig$values < 0] <- 0
    
    # Reconstituindo a matriz corrigida
    B_corrigida <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    
    return(B_corrigida)
  }
  
  # Corrigir matriz B
  b_corrigida <- corrigir_com_pca(b_manova)
  
  print(eigen(b_corrigida)$values)
  
  # Inspecionar a matriz B
  print(b_manova)
  
  cov_para_cor(b_corrigida)
  
  # Checar os autovalores para avaliar estabilidade
  autovalores <- eigen(b_manova)$values
  print(autovalores)
  
  corrigir_autovalores <- function(B) {
    eig <- eigen(B)
    eig$values[eig$values < 0] <- 0
    return(eig$vectors %*% diag(eig$values) %*% t(eig$vectors))
  }
  
  # Corrigir matriz B
  b_corrigida <- corrigir_autovalores(b_manova)
  eigen(b_corrigida)$values
  
  qr(b_corrigida)$rank
  
  # Gerar matriz de correlação
  R_b <- cov_para_cor(b_manova)
  print(R_b)
  
  
  
  medias = rbind(colMeans(sub_sexes_f[,50:88]), colMeans(sub_sexes_m[,50:88]))

  bes_medias[[i]] = simetrizar_matriz(cov(medias))
  
  names(bes)[[i]] = paste(genus, sp)
  names(bes_medias)[[i]] = paste(genus, sp)
}

save(bes, file = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Matriz_B.RData")
save(bes_medias, file = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Matriz_B_medias.RData")

#
load("~/Dropbox/Doc/Code/evowm/R/Novo_Output/Matriz_B.RData")
load("~/Dropbox/Doc/Code/evowm/R/Novo_Output/Matriz_B_medias.RData")
data = readRDS("~/Dropbox/Doc/Code/evowm/R/Novo_Output/Averages_PCS_Extant_Species.RData")

# Função para converter matriz de covariância em matriz de correlação


mandril = 








# correlações maximas?
RandomSkewers(bes[[1]], bes_medias[[1]])

#
x = eigen(bes_medias[[2]])$values
y = eigen(bes[[2]])$values

x == y

#
is.positive.definite(bes[[30]])

for(i in 1:length(bes)){
  check = any(eigen(bes[[i]])$values < 0)
  print(check)
}

solve(as.matrix(bes[[1]]))

#
gibao = bes$`Hylobates muelleri` / ((data$NoLog_Averages$Hylobates_muelleri$Machos + data$NoLog_Averages$Hylobates_muelleri$Fêmeas) / 2)

D = diag(1/sqrt(diag(gibao)))
gibao_2 = D%*%gibao%*%D

mandril = bes$`Mandrillus sphinx` / ((data$NoLog_Averages$Mandrillus_sphinx$Machos + data$NoLog_Averages$Mandrillus_sphinx$Fêmeas) / 2)

D = diag(1/sqrt(diag(mandril)))
mandril_2 = D%*%mandril%*%D

gib <- matrix(as.numeric(gibao_2), nrow = nrow(gibao_2))
is.numeric(gib)
sd(gib)
cv = (sd(gib) / mean(gib))
cat("Coeficiente de Variação (CV):", round(cv, 2), "%\n")


mand <- matrix(as.numeric(mandril), nrow = nrow(mandril))
is.numeric(mand)
sd(mand)
cv = (sd(mand) / mean(mand))
cat("Coeficiente de Variação (CV):", round(cv, 2), "%\n")


