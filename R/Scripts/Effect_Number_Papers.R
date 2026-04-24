# number of papers effect on the variation of behaviours
# testing if more papers leads to more 'vaiation' in classes

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts")

library(stringr)
library(MASS)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  sum(rp^2) / rdf
}

# read mating data
dados = read.csv("~/Dropbox/Doc/Data/wos_mating_systems/Matings.csv")
dados = dados[!is.na(dados$DOMINANCE) & !is.na(dados$AGGRESSION), ]

# get species names
especies = vector()
for(i in 1:nrow(dados)){
  especies[i] = paste(dados$GENUS[i], dados$SPECIES[i], sep = "_")
}

dados$n_classes_social <- str_count(dados$SOCIAL_ORGANIZATION, "/") + 1
dados$n_classes_matings <- str_count(dados$MATING_SYSTEM, "/") + 1
dados$n_classes_dominance <- str_count(dados$DOMINANCE, "/") + 1
dados$n_classes_aggression <- str_count(dados$AGGRESSION, "/") + 1

m_social  <- glm(n_classes_social ~ N_PAPERS,
                 family = poisson(link = "log"),
                 data = dados)

m_mating  <- glm(n_classes_matings ~ N_PAPERS,
                 family = poisson(link = "log"),
                 data = dados)

m_dominance  <- glm(n_classes_dominance ~ N_PAPERS,
                 family = poisson(link = "log"),
                 data = dados)

m_aggression  <- glm(n_classes_aggression ~ N_PAPERS,
                    family = poisson(link = "log"),
                    data = dados)

lm(PROP_MALES_FEMALES ~ N_PAPERS, data = dados)
summary(lm(PROP_MALES_FEMALES ~ N_PAPERS, data = dados))

summary(m_social)
summary(m_mating)
summary(m_dominance)
summary(m_aggression)

overdisp_fun(m_social)
overdisp_fun(m_mating)
overdisp_fun(m_dominance)
overdisp_fun(m_aggression)

