#CorVec_Dimorphism_PCs
# vector correlation between sexual dimorphism and 
# Principal Components PCs for all species

# load packages and functions
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}

geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}


prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# load data to correlation
medidas = readRDS("~/Dropbox/Novo_Output/Averages_PCS_Extant_Species.RData")

# read all VCV matrices
# read all vcv matrices
setwd("~/Dropbox/Doc/Output/nolog_vcv")
temp = list.files(pattern = "*.txt")

# remove species which are not in the phylogeny
to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

index = grep(paste(to_remove, collapse = "|"), temp)
temp = temp[-index]

vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
names(vcv)  = gsub(".txt", replacement= "", temp)

# get normalized dimorphism
dimor = vector()
for(i in 1:length(medidas$NoLog_Averages)){
  sp = medidas$NoLog_Averages[[i]]
  dimorf = unlist(sp[1]) - unlist(sp[2])
  
  tamanho_macho = medidas$NoLog_ByTrait_Averages[[i]]$Machos
  tam_cra = geomean(tamanho_macho)
  
  dimor[i] = dimorf / tam_cra
}

# loop
align = vector()
for(i in 1:length(medidas$NoLog_Averages)){
  sp = medidas$Species[[i]]
  med = medidas$NoLog_ByTrait_Averages[[sp]]
  covar = vcv[[sp]]
  
  # pc
  pc_1 = eigen(covar)$vectors[,1]
  
  # dimorphism
  d = (medidas$NoLog_ByTrait_Averages[[sp]]$Machos - medidas$NoLog_ByTrait_Averages[[sp]]$Fêmeas) / 
    geomean(medidas$NoLog_ByTrait_Averages[[i]]$Machos)
  
  align[i] = abs(corVector(pc_1, d))
  
}

# read evolvability
evolva = readRDS("~/Dropbox/Novo_Output/Evolvability.RData")

toplot = data.frame(medidas$Species, unlist(evolva$Standart_Evolvability), align)

#saveRDS(toplot, file = "~/Dropbox/Novo_Output/CorVec_Dimorphism_PCs.RData")

p = ggplot(data = toplot) +
  geom_point(aes(x = align, y = unlist.evolva.Standart_Evolvability.)) +
  geom_smooth(aes(x = align, y = unlist.evolva.Standart_Evolvability.), method = "lm", se = FALSE, formula = y ~ x) +
  xlab("Alignment between dimorphism and Pmax") +
  ylab("Standardized Evolvability") +
  theme_classic()

p

#cor(log(toplot$align), log(toplot$unlist.evolva.Standart_Evolvability.))

ggsave(p, filename = "~/Dropbox/Novo_Output/Alignment_Evolvability.pdf", dpi = 600,
       width = 35, height = 20, units = "cm",  bg = "transparent")
