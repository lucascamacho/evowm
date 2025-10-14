#
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(MASS)){install.packages("MASS"); library(MASS)}
if(!require(GGally)){install.packages("GGally"); library(GGally)}

medidas = read.table("~/Dropbox/cor_PCS_dimorphism_extant.csv", header = TRUE, sep = ",")
matings = readRDS("~/Dropbox/New_Haplorrhini_MDS_Matings.RDS") 

# Primeiro, armazenamos a ordem que precisamos
ordem = match(matings$especies, medidas$matings.especies)

# Reorganiza as linhas de medidas_simi_filtrado
medidas = medidas[ordem, ]

############################
# Haplorrhini
data = data.frame(medidas, matings)

vars <- c("dimor", "dados.SOCIAL_ORGANIZATION", "dados.MATING_SYSTEM", "dados.PROP_MALES_FEMALES", 
          "dados.DOMINANCE", "dados.AGGRESSION", "dados.N_PAPERS")


ggpairs(data[, vars])

# Arrumar align para testes
align_z = 0.5 * log((1 + data$align_1) / (1 - data$align_1))
data = cbind(data, align_z)

# plot dimorphism x align PC1 Dimor separando os grupos de comportamento. Escolher qual variável.
ggplot(data, aes(x = align_z, y = dimor, color = as.factor(dados.MATING_SYSTEM))) +
  geom_point(size = 2) +
  stat_ellipse(type = "norm", level = 0.95) +
  theme_minimal() +
  labs(color = "Grupo", x = "Align PC1-Dimor (z de Fisher)", y = "Dimorfismo",
       title = "Distribuição dos grupos com elipses de confiança")

# multivariate test
modelo = manova(cbind(Dimorphism, align_z) ~ scores.fit....1. + scores.fit....2. +
                  scores.fit....1.:scores.fit....2., data = data)

summary(modelo)       # Testes multivariados
summary.aov(modelo)   # ANOVAs univariadas

###################
# Platyrrhini
data = data.frame(medidas_simi_filtrado, matings_filtrado)
index = which(data$dados.PARVORDER == "Platyrrhini")
data = data[index,]

vars <- c("Dimorphism", "dados.SOCIAL_ORGANIZATION", "dados.MATING_SYSTEM", "dados.PROP_MALES_FEMALES", 
          "dados.DOMINANCE", "dados.AGGRESSION")

ggpairs(data[, vars])

# Arrumar align para testes
align_z = 0.5 * log((1 + data$cor_PC1_Dimor) / (1 - data$cor_PC1_Dimor))
data = cbind(data, align_z)

# plot dimorphism x align PC1 Dimor separando os grupos de comportamento. Escolher qual variável.
ggplot(data, aes(x = align_z, y = Dimorphism, color = as.factor(dados.SOCIAL_ORGANIZATION))) +
  geom_point(size = 2) +
  stat_ellipse(type = "norm", level = 0.95) +
  theme_minimal() +
  labs(color = "Grupo", x = "Align PC1-Dimor (z de Fisher)", y = "Dimorfismo",
       title = "Distribuição dos grupos com elipses de confiança")

# multivariate test
modelo = manova(cbind(Dimorphism, align_z) ~ scores.fit....1. + scores.fit....2. +
                  scores.fit....1.:scores.fit....2., data = data)

summary(modelo)       # Testes multivariados
summary.aov(modelo)   # ANOVAs univariadas

#############
# Catarrhini
data = data.frame(medidas_simi_filtrado, matings_filtrado)
index = which(data$dados.PARVORDER == "Catarrhini")
data = data[index,]

vars <- c("Dimorphism", "dados.SOCIAL_ORGANIZATION", "dados.MATING_SYSTEM", "dados.PROP_MALES_FEMALES", 
          "dados.DOMINANCE", "dados.AGGRESSION")

ggpairs(data[, vars])

# Arrumar align para testes
align_z = 0.5 * log((1 + data$cor_PC1_Dimor) / (1 - data$cor_PC1_Dimor))
data = cbind(data, align_z)

# plot dimorphism x align PC1 Dimor separando os grupos de comportamento. Escolher qual variável.
ggplot(data, aes(x = align_z, y = Dimorphism, color = as.factor(dados.SOCIAL_ORGANIZATION))) +
  geom_point(size = 2) +
  stat_ellipse(type = "norm", level = 0.95) +
  theme_minimal() +
  labs(color = "Grupo", x = "Align PC1-Dimor (z de Fisher)", y = "Dimorfismo",
       title = "Distribuição dos grupos com elipses de confiança")

# multivariate test
modelo = manova(cbind(Dimorphism, align_z) ~ scores.fit....1. + scores.fit....2. +
                  scores.fit....1.:scores.fit....2., data = data)

summary(modelo)       # Testes multivariados
summary.aov(modelo)   # ANOVAs univariadas
