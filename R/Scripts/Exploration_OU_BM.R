setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(mvSLOUCH)
library(ape)
library(dplyr)
library(phytools)
library(phytools)
library(gt)

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

#Load and plot the phylogeny
# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
species = mds$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

align_filtrado <- align[align$matings.especies %in% mds$especies, ]
evolvas_filtrado <- evolvas[evolvas$species %in% mds$especies, ]

# RESULTS AICC
# Catarrhini
cata = load("~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

# See AICc values
aicc_catar = c(fits$BM$ParamSummary$aic.c,
               fits$OU_SIMPLE$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL$FinalFound$ParamSummary$aic.c,
               fits$OU_MATINGS$FinalFound$ParamSummary$aic.c,
               fits$OU_AGGRESSION$FinalFound$ParamSummary$aic.c,
               fits$OU_INTEGRATION$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL_MATINGS$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL_MATINGS_AGGRESSION$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL_MATINGS_AGGRESSION_INTEGRATION$FinalFound$ParamSummary$aic.c)

names(aicc_catar) = c("BM",
                      "SIMPLE_OU",
                      "OU_SOCIAL",
                      "OU_MATINGS",
                      "OU_AGGRESSION",
                      "OU_INTEGRATION",
                      "OU_SOCIAL_MATINGS",
                      "SOCIAL_MATINGS_AGGRESSION",
                      "SOCIAL_MATINGS_AGGRESSION_INTEGRATION")

# -----------------------------
# 2) Função para montar tabela
# -----------------------------
make_aicc_table <- function(aicc_values) {
  data.frame(
    model = names(aicc_values),
    AICc = as.numeric(aicc_values)
  ) %>%
    arrange(AICc) %>%
    mutate(
      deltaAICc = AICc - min(AICc, na.rm = TRUE),
      weight = exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc))
    )
}

tab_catar <- make_aicc_table(aicc_catar)

cat("\n=== Catarrhini ===\n")
print(tab_catar, digits = 3)

######################
# Platyrrini
platy = load("~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fits.RData")

# See AICc values
aicc_platy = c(fits$BM$ParamSummary$aic.c,
               fits$OU_SIMPLE$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL$FinalFound$ParamSummary$aic.c,
               fits$OU_MATINGS$FinalFound$ParamSummary$aic.c,
               fits$OU_PROP$FinalFound$ParamSummary$aic.c,
               fits$OU_AGGRESSION$FinalFound$ParamSummary$aic.c,
               fits$OU_INTEGRATION$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL_MATINGS$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL_MATINGS_PROP$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL_MATINGS_PROP_AGGRESSION$FinalFound$ParamSummary$aic.c,
               fits$OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION$FinalFound$ParamSummary$aic.c
               )

names(aicc_platy) = c("BM",
                      "SIMPLE_OU",
                      "OU_SOCIAL",
                      "OU_MATINGS",
                      "OU_PROP",
                      "OU_AGGRESSION",
                      "OU_INTEGRATION",
                      "OU_SOCIAL_MATINGS",
                      "OU_SOCIAL_MATINGS_PROP",
                      "OU_SOCIAL_MATINGS_PROP_AGGRESSION",
                      "OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION"
                      )


# -----------------------------
# 2) Função para montar tabela
# -----------------------------
make_aicc_table <- function(aicc_values) {
  data.frame(
    model = names(aicc_values),
    AICc = as.numeric(aicc_values)
  ) %>%
    arrange(AICc) %>%
    mutate(
      deltaAICc = AICc - min(AICc, na.rm = TRUE),
      weight = exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc))
    )
}

tab_platy <- make_aicc_table(aicc_platy)

cat("\n=== Platyrrhini ===\n")
print(tab_platy, digits = 3)

###### Half lifes
cata = load("~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

hist(fits$OU_SIMPLE$FinalFound$ParamSummary$phyl.halflife$halflive[2,1:39])

hist(diag(fits$OU_SIMPLE$FinalFound$ParamSummary$expmtA))


platy = load("~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fits.RData")

hist(fits$OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION$FinalFound$ParamSummary$phyl.halflife$halflive[2,1:39])

hist(diag(fits$OU_SOCIAL_MATINGS_PROP_AGGRESSION_INTEGRATION$FinalFound$ParamSummary$expmtA))
