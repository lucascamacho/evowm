setwd("~/Dropbox/Doc/Doc/Code/evowm/R/Scripts/")

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

# RESULTS
# Catarrhini
cata = load("~/Dropbox/catarrhini_mvslouch_fits.RData")

catar <- fits$OU_MDS2$FinalFound$ParamsInModel
catar_df <- tibble(
  Clade = "Catarrhini",
  alpha = catar$A[1,1],
  #theta_NMDS1 = catar$B[1,1],
  theta_NMDS1 = NULL,
  theta_NMDS2 = catar$B[1,1],
  sigma2 = catar$Syy[1,1],
  mPsi = catar$mPsi[1,1],
  vY0 = catar$vY0[1,1],
  beta1 = catar$B[1,1],
  #beta2 = catar$B[1,2],
  beta2 = NULL,
  intercept = as.numeric(catar$mPsi),
)


# See AICc values
aicc_catar = c(fits$BM$ParamSummary$aic.c,
               fits$Simple_OU$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS1$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS2$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS12$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS12_size$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS12_size_cond_evolvability$FinalFound$ParamSummary$aic.c,
               fits$OU_full$FinalFound$ParamSummary$aic.c)
names(aicc_catar) = c("BM",
                      "Simple_OU",
                      "OU_MDS1",
                      "OU_MDS2",
                      "OU_MDS12",
                      "OU_MDS12_size",
                      "OU_MDS12_size_cond_evolvability",
                      "OU_full")

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

# Platyrrini
platy = load("~/Dropbox/platyrrhini_mvslouch_fits.RData")

platy <- fits$OU_MDS12$FinalFound$ParamsInModel

platy_df <- tibble(
  Clade = "Platyrrhini",
  alpha = platy$A[1,1],
  theta_NMDS1 = platy$B[1,1],
  theta_NMDS2 = platy$B[1,2],
  sigma2 = platy$Syy[1,1],
  mPsi = platy$mPsi[1,1],
  vY0 = platy$vY0[1,1], 
  beta1 = platy$B[1,1],
  beta2 = platy$B[1,2],
  intercept = as.numeric(platy$mPsi),)

# See AICc values
aicc_platy = c(fits$BM$ParamSummary$aic.c,
               fits$Simple_OU$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS1$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS2$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS12$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS12_size$FinalFound$ParamSummary$aic.c,
               fits$OU_MDS12_size_cond_evolvability$FinalFound$ParamSummary$aic.c,
               fits$OU_full$FinalFound$ParamSummary$aic.c)
names(aicc_platy) = c("BM",
                      "Simple_OU",
                      "OU_MDS1",
                      "OU_MDS2",
                      "OU_MDS12",
                      "OU_MDS12_size",
                      "OU_MDS12_size_cond_evolvability",
                      "OU_full")

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


### 
params_df <- bind_rows(platy_df, catar_df)


# -----------------------------
# Criar tabela gt
# -----------------------------
params_df %>%
  mutate(
    half_life = ifelse(alpha > 0, log(2)/alpha, Inf),
    n_species = c(38, 24)  # primeira linha Catarrhini, segunda Platyrrhini
  ) %>%
  gt() %>%
  tab_header(
    title = "Parameters Model OU - Catarrhini e Platyrrhini"
  ) %>%
  fmt_number(
    columns = tidyselect::all_of(c("alpha", "theta_NMDS1", "theta_NMDS2",
                                   "sigma2", "mPsi", "vY0", "half_life")),
    decimals = 4
  ) %>%
  cols_label(
    alpha = "α (Alpha)",
    theta_NMDS1 = "θ NMDS1",
    theta_NMDS2 = "θ NMDS2",
    sigma2 = "σ² (Syy)",
    mPsi = "mΨ",
    vY0 = "vY0",
    half_life = "t½ (Half-Life)",
    n_species = "N Species",
    beta1 = "β NMDS1",
    beta2 = "β NMDS2"
  )
