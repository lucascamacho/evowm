# Evolvability of sexual dimorphism
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(ggplot2)
library(openxlsx)
library(evolqg)
library(evolvability)

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# read V/CV matrices
setwd("~/Dropbox/Doc/Data/vcv/")
temp = list.files(pattern="*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement= "", temp)

#
matings = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

results <- list()  # lista vazia
for(i in 1:length(matings$especies)){
  #
  sp = matings$especies[i]
  medidas = medias$ByTrait_Averages[[sp]]
  
  # se não tem medidas, pula
  if(is.null(medidas)){
    next
  }
  
  # se não está no vcv, pula
  if(!(sp %in% names(vcv))){
    next
  }
  
  # cálculos
  dimor <- medidas$Machos - medidas$Fêmeas
  tam_cra <- geomean(medidas$Machos)
  size <- (geomean(medidas$Machos) + geomean(medidas$Fêmeas)) / 2
  
  # calcula médias dos traits
  medias_por_trait <- (medidas$Machos + medidas$Fêmeas) / 2
  medias_por_trait <- as.numeric(medias_por_trait)   # garante vetor numérico
  
  # padronizar matriz
  covar <- vcv[[sp]]
  covar <- as.matrix(covar)                          # garante matriz numérica
  mean_prod <- outer(medias_por_trait, medias_por_trait, "*")   # matriz com produtos das médias
  mat_scaled <- covar / mean_prod   # padroniza tudo
  diag(mat_scaled) <- diag(covar) / medias_por_trait   # substitui a diagonal (variâncias / média)
  
  # normalizações
  dimor = dimor / size
  dimor_norm <- dimor / sqrt(sum(dimor ^ 2))
  dimorfism <- mean(dimor)
  
  # evolvability
  resposta <- as.matrix(covar) %*% as.vector(dimor_norm)
  evolv <- sum(resposta * dimor_norm) / (size * size)
  
  # std evolvability
  std_evolv <- sqrt(evolv) / (size * size)
  
  # integration OK
  rownames(covar) <- colnames(covar) <- NULL
  integ <- CalcEigenVar(as.matrix(covar))
  
  # avg_evolvability OK
  avg_evolvability = sum(diag(as.matrix(covar))) / qr(as.matrix(covar))$rank
  
  # salva numa linha
  results[[sp]] <- data.frame(
    species = sp,
    Evolvability = evolv,
    Standart_Evolvability = std_evolv,
    Dimorfism = dimorfism,
    Size = size,
    Integration = integ,
    Average_Evolvability = avg_evolvability,
    Conditional_Evolvability = (evolvabilityBeta(as.matrix(covar), as.vector(dimor_norm))$c) / (size * size),
    Integration_A = 1 - evolvabilityBeta(as.matrix(covar), as.vector(dimor_norm))$a
  )
}

# junta tudo num único data.frame
evolvas <- do.call(rbind, results)
rownames(evolvas) <- NULL

# save results
saveRDS(evolvas, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")

evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")

# plot
evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Size, y = Standart_Evolvability)) +
  xlab("Size") +
  ylab("SD_Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot


evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Average_Evolvability, y = Evolvability)) +
  xlab("Avg_Evolvability") +
  ylab("SD_Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Evolvability, y = Dimorfism)) +
  xlab("Evolvability") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot


evolv_size_plot = ggplot(data = evolvas) +
  #geom_point(aes(x = Size, y = Evolvability)) +
  geom_text(aes(x = Size, y = Evolvability, label = species)) +
  xlab("Size") +
  ylab("Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_size_plot

integ_sd_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Integration, y = Dimorfism)) +
  #  geom_hline(yintercept = mean(as.numeric(std_evolvability)), linetype = "dashed") +
  xlab("Integration_A") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

integ_sd_plot

integ_evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Integration, y = Evolvability)) +
  #  geom_hline(yintercept = mean(as.numeric(std_evolvability)), linetype = "dashed") +
  xlab("Integration_A") +
  ylab("Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

integ_evolv_plot

evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Evolvability, y = Dimorfism, label = species)) +
  xlab("Evolvability") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Conditional_Evolvability, y = Dimorfism, label = species)) +
  xlab("Conditional_Evolvability") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Conditional_Evolvability, y = Evolvability, label = species)) +
  xlab("Conditional_Evolvability") +
  ylab("Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Conditional_Evolvability, y = Integration_A, label = species)) +
  xlab("Conditional_Evolvability") +
  ylab("Integration_A") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

#
alinhamentos = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv")

dimor_data = data.frame(evolvas$species, evolvas$Standart_Evolvability, evolvas$Dimorfism, 
                        evolvas$Conditional_Evolvability, evolvas$Integration_A,
                        alinhamentos$align_1)

dimor_plot = ggplot(data = dimor_data) +
  geom_point(aes(x = evolvas.Standart_Evolvability, y = evolvas.Dimorfism)) +
  geom_smooth(aes(x = evolvas.Standart_Evolvability, y = evolvas.Dimorfism), method = "lm", se = FALSE, formula = y ~ x) +
  xlab("Standardized Evolvability") +
  ylab("Standardized Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),  # Tamanho das legendas dos eixos
    axis.text = element_text(size = 12)    # Tamanho dos números dos eixos
  )

dimor_plot

#ggsave(dimor_plot, filename = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Dimorphism_Evolvability.pdf", dpi = 600,
#       width = 35, height = 20, units = "cm",  bg = "transparent")

dimor_plot_2 = ggplot(data = dimor_data) +
  geom_point(aes(x = alinhamentos.align_1, y = evolvas.Dimorfism)) +
  geom_smooth(aes(x = alinhamentos.align_1, y = evolvas.Dimorfism), method = "lm", se = FALSE, formula = y ~ x) +
  xlab("Pmax x Sexual Dimorphism Vector") +
  ylab("Standardized Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),  # Tamanho das legendas dos eixos
    axis.text = element_text(size = 12)    # Tamanho dos números dos eixos
  )

dimor_plot_2

#ggsave(dimor_plot_2, filename = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Alinhamento_Dimorphism_Evolvability.pdf", dpi = 600,
#       width = 35, height = 20, units = "cm",  bg = "transparent")
