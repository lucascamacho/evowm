# Evolvability of sexual dimorphism

library(ggplot2)
library(openxlsx)

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

#read V/CV matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Novo_Output/p_vcv_gabriel/catarrhini/")
temp = list.files(pattern="*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement= "", temp)


# read all vcv matrices
#setwd("~/Dropbox/Doc/Output/nolog_vcv")
#temp = list.files(pattern = "*.txt")

# remove species which are not in the phylogeny
#to_remove = c("Bunopithecus_hoolock", "Cercopithecus_lhoesti", "Cercopithecus_preussi",
#              "Kasi_johnii", "Kasi_vetulus", "Lophocebus_opdenboschi")

#index = grep(paste(to_remove, collapse = "|"), temp)
#temp = temp[-index]

#vcv = lapply(temp, read.table, header = TRUE, row.names = 1)
#names(vcv)  = gsub(".txt", replacement= "", temp)

# read datas
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Novo_Output/Averages_PCS_Extant_Species.RData")

# create lists for results
evolvability = list()
std_evolvability = list()
c_evolvability = list()
std_c_evolvability = list()
for(i in 1:length(names(vcv))){
  #choose species, measures and VCV
  sp = names(vcv)[i]
  medidas = medias$NoLog_ByTrait_Averages[[sp]]
  
  if(is.null(medidas) == TRUE){
    next
  }
    
  
  dimor = medidas$Machos - medidas$Fêmeas
  
  tam_cra = geomean(medidas$Machos)
  
  covar = vcv[[sp]]
  
  # evolvability
  dimor_norm = dimor / sqrt(sum(dimor ^ 2))

  resposta = as.matrix(covar) %*% as.vector(dimor_norm)

  evolv = sum(resposta * dimor_norm)
  evolvability[[i]] = evolv
  
  # std evolvability
  evolv_norm = sqrt(evolv) / tam_cra
  std_evolvability[[i]] = evolv_norm

}

# rename list of results
evolvas = list(evolvability, std_evolvability)

names(evolvas) = c("Evolvability", "Standart_Evolvability")

names(evolvas$Evolvability) = names(vcv)
names(evolvas$Standart_Evolvability) = names(vcv)

# save results
#saveRDS(evolvas, file = "~/Dropbox/Novo_Output/Evolvability.RData")

# prepare data frame
toplot = data.frame(names(vcv)[-8], unlist(evolvas$Evolvability), unlist(evolvas$Standart_Evolvability))
names(toplot) = c("species", "evolvability", "std_evolvability")

# plot
evolv_plot = ggplot(data = toplot) +
  geom_point(aes(x = species, y = as.numeric(std_evolvability))) +
  geom_hline(yintercept = mean(as.numeric(std_evolvability)), linetype = "dashed") +
  xlab("Species") +
  ylab("Standardized Evolvability") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

evolv_plot

# save the final plot
#ggsave(evolv_plot, filename = "~/Dropbox/Novo_Output/Evolvability.pdf", dpi = 600,
#       width = 65, height = 35, units = "cm",  bg = "transparent")

dimorfismos = vector()
alinhamentos = vector()
for(i in 1:length(names(vcv))){
  sp = names(vcv)[i]
  sp = medias$NoLog_ByTrait_Averages[[sp]]
  
  if(is.null(sp) == TRUE){
    next
  }
  
  dimorf = unlist(sp[1]) - unlist(sp[2])
  
  tamanho_macho = unlist(sp[1])
  tam_cra = geomean(tamanho_macho)
    
  dimorfismos[i] = mean(dimorf) / tam_cra
  
  #
  sp = names(vcv)[i]
  covar = vcv[[sp]]
  
  p_max = eigen(covar)$vectors[,1]
  alinhamentos[i] = abs(corVector(p_max, dimorf))
}

dimor_data = data.frame(toplot$species, toplot$std_evolvability, na.omit(dimorfismos), na.omit(alinhamentos))

dimor_plot = ggplot(data = dimor_data) +
  geom_point(aes(x = toplot.std_evolvability, y = na.omit.dimorfismos.)) +
  geom_smooth(aes(x = toplot.std_evolvability, y = na.omit.dimorfismos.), method = "lm", se = FALSE, formula = y ~ x) +
  xlab("Standardized Evolvability") +
  ylab("Standardized Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),  # Tamanho das legendas dos eixos
    axis.text = element_text(size = 12)    # Tamanho dos números dos eixos
  )

dimor_plot

ggsave(dimor_plot, filename = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Dimorphism_Evolvability.pdf", dpi = 600,
       width = 35, height = 20, units = "cm",  bg = "transparent")

dimor_plot_2 = ggplot(data = dimor_data) +
  geom_point(aes(x = na.omit.alinhamentos., y = na.omit.dimorfismos.)) +
  geom_smooth(aes(x = na.omit.alinhamentos., y = na.omit.dimorfismos.), method = "lm", se = FALSE, formula = y ~ x) +
  xlab("Pmax x Sexual Dimorphism Vector") +
  ylab("Standardized Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),  # Tamanho das legendas dos eixos
    axis.text = element_text(size = 12)    # Tamanho dos números dos eixos
  )

dimor_plot_2

cor(dimor_data$na.omit.dimorfismos., dimor_data$na.omit.alinhamentos.)

ggsave(dimor_plot_2, filename = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Alinhamento_Dimorphism_Evolvability.pdf", dpi = 600,
       width = 35, height = 20, units = "cm",  bg = "transparent")
