# PGLS POR ESPECIE
# CATARRHINI, PLATYRRHINI, HAPLORRHINI
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(slouch)
library(ape)
library(dplyr)
library(caper)
library(phytools)
library(ggplot2)
library(ggrepel)
library(scales) 

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/ancestrals_averages_PCS_autovalues_primates.RDS")
#nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
  #nicho$`esp3$...1`,
  align$matings.especies,
  evolvas$Species,
  mds$especies
))

#nicho_f   <- nicho[nicho$`esp3$...1` %in% sp_comuns, ]
align_f   <- align[align$matings.especies %in% sp_comuns, ]
evolvas_f <- evolvas[evolvas$Species %in% sp_comuns, ]
mds_f     <- mds[mds$especies %in% sp_comuns, ]

# read all VCV matrices
setwd("~/Dropbox/Doc/Data/vcv/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement = "", temp)
vcv$Lagothrix_lagothricha <- as.data.frame(lapply(vcv$Lagothrix_lagothricha, function(x) as.numeric(trimws(x))))
vcv$Cacajao_calvus = read.csv("~/Dropbox/Doc/Data/p_vcv_gabriel/Cacajao_calvus.csv", header = FALSE, sep = ";", dec = ",")
vcv$Cacajao_calvus <- as.data.frame(lapply(vcv$Cacajao_calvus, function(x) as.numeric(trimws(x))))

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
tree$tip.label[tree$tip.label == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"
species = mds_f$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

#
# media_especie <- nicho %>%
#   dplyr::group_by(.data$`esp3$...1`) %>%
#   dplyr::summarise(
#   climacv_mean = mean(.data$climacv, na.rm = TRUE),
#   CVn_mean = mean(.data$CVn, na.rm = TRUE)
#   )

#nicho_f <- media_especie[match(tree$tip.label, media_especie$`esp3$...1`), ]
evolvas_f <- evolvas_f[match(tree$tip.label, evolvas_f$Species), ]
align_f <- align_f[match(tree$tip.label, align_f$matings.especies), ]
mds_f <- mds_f[match(tree$tip.label, mds_f$especies), ]

################################################################################
# A PARTIR DAQUI
# Catarrhini #1:38 clima e cvn, 1:34 com clima e cvn
align_c = align_f[1:38,]
evolvas_c = evolvas_f[1:38,]
mds_c = mds_f[1:38,]
tree_c = drop.tip(tree, setdiff(tree$tip.label, align_c$matings.especies))
#nicho_c = nicho_f[1:34,]

#
df <- data.frame(
  species = mds_c$especies,
  parvorder = mds_c$dados.PARVORDER,
  social = mds_c$dados.SOCIAL_ORGANIZATION,
  mating = mds_c$dados.MATING_SYSTEM,
  prop = mds_c$dados.PROP_MALES_FEMALES,
  dominance = mds_c$dados.DOMINANCE,
  aggression = mds_c$dados.AGGRESSION,
  dimorfism = evolvas_c$Dimorfism,
  integration = evolvas_c$Integration,
  size = evolvas_c$Size,
  nmds1 = mds_c$scores.fit....1.,
  nmds2 = mds_c$scores.fit....2.,
  normas = align_c$normas
  #clima = nicho_c$climacv_mean,
  #cvn = nicho_c$CVn_mean
)

df <- df[match(tree_c$tip.label, df$species), ]

predictors <- c("integration", "size", "nmds1", "nmds2")
response <- "normas"

rownames(df) = df$species  # ou use df$genus se sua árvore for por gênero

# Cria o objeto comparative.data
comp <- comparative.data(
  phy = tree_c,       # sua árvore
  data = df,          # seu dataframe
  names.col = "species", # coluna com os nomes que batem com a árvore
  vcv = TRUE,         # para calcular a matriz de covariância
  na.omit = FALSE     # mantém NAs se houver (ou TRUE para remover)
)


# Lista vazia para fórmulas
formulas <- list()

# Loop para gerar combinações de 1 a 4 preditores
for (n in 1:length(predictors)) {
  combs <- combn(predictors, n, simplify = FALSE)
  for (x in combs) {
    f <- paste(response, "~", paste(x, collapse = " + "))
    formulas <- c(formulas, f)
  }
}

# pgls
results <- lapply(formulas, function(f){
  fit <- pgls(as.formula(f), data = comp)
  aic = AIC(fit)
  # número de observações
  n <- length(fit$residuals)
  # número de parâmetros (intercepto + preditores)
  k <- length(coef(fit))
  # AICc manual
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)  
  r2 <- summary(fit)$adj.r.squared
  coef <- summary(fit)$coefficients
  list(formula = f, fit = fit, AICc = aicc, R2 = r2, coefficients = coef)
})

results_df <- bind_rows(lapply(results, function(x){
  data.frame(
    formula = x$formula,
    AICc = x$AICc,
    R2 = x$R2
  )
}))

# Ordenar por AICc
results_df$deltaAICc = results_df$AICc - min(results_df$AICc)
results_df$weight = exp(-0.5 * results_df$deltaAICc)
results_df$weight = results_df$weight / sum(results_df$weight)
results_df <- results_df %>% arrange(deltaAICc)
head(results_df)

#write.csv(results_df, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_PGLS.csv")

best_cat <- pgls(normas ~ size, data = comp)
summary(best_cat)

plot(fitted(best_cat), residuals(best_cat),
     xlab = "Valores ajustados",
     ylab = "Resíduos",
     main = "Catarrhini: resíduos vs ajustados")
abline(h = 0, col = "red", lty = 2)

hist(residuals(best_cat), main="Catarrhini: Histograma de resíduos", xlab="Resíduos")
qqnorm(residuals(best_cat)); qqline(residuals(best_cat), col="red")

shapiro.test(residuals(best_cat))

res_cat <- residuals(best_cat)
lambda_cat <- phylosig(tree_c, res_cat, method = "lambda")
lambda_cat

#################################################################################
# Platyrrhini # 39:65 sem clima e cvn, 35:53 com clima e cvn
align_p = align_f[39:65,]
evolvas_p = evolvas_f[39:65,]
mds_p = mds_f[39:65,]
tree_p = drop.tip(tree, setdiff(tree$tip.label, align_p$matings.especies))
#nicho_p = nicho_f[35:53,]

#
df <- data.frame(
  species = mds_p$especies,
  parvorder = mds_p$dados.PARVORDER,
  social = mds_p$dados.SOCIAL_ORGANIZATION,
  mating = mds_p$dados.MATING_SYSTEM,
  prop = mds_p$dados.PROP_MALES_FEMALES,
  dominance = mds_p$dados.DOMINANCE,
  aggression = mds_p$dados.AGGRESSION,
  dimorfism = evolvas_p$Dimorfism,
  integration = evolvas_p$Integration,
  size = evolvas_p$Size,
  nmds1 = mds_p$scores.fit....1.,
  nmds2 = mds_p$scores.fit....2.,
  normas = align_p$normas
  #clima = nicho_p$climacv_mean,
  #cvn = nicho_p$CVn_mean
)

df <- df[match(tree_p$tip.label, df$species), ]

predictors <- c("integration", "size", "nmds1", "nmds2")
response <- "normas"

rownames(df) = df$species  # ou use df$genus se sua árvore for por gênero

# Cria o objeto comparative.data
comp <- comparative.data(
  phy = tree_p,       # sua árvore
  data = df,          # seu dataframe
  names.col = "species", # coluna com os nomes que batem com a árvore
  vcv = TRUE,         # para calcular a matriz de covariância
  na.omit = FALSE     # mantém NAs se houver (ou TRUE para remover)
)


# Lista vazia para fórmulas
formulas <- list()

# Loop para gerar combinações de 1 a 4 preditores
for (n in 1:length(predictors)) {
  combs <- combn(predictors, n, simplify = FALSE)
  for (x in combs) {
    f <- paste(response, "~", paste(x, collapse = " + "))
    formulas <- c(formulas, f)
  }
}

# pgls
results <- lapply(formulas, function(f){
  fit <- pgls(as.formula(f), data = comp)
  aic = AIC(fit)
  # número de observações
  n <- length(fit$residuals)
  # número de parâmetros (intercepto + preditores)
  k <- length(coef(fit))
  # AICc manual
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)  
  r2 <- summary(fit)$adj.r.squared
  coef <- summary(fit)$coefficients
  list(formula = f, fit = fit, AICc = aicc, R2 = r2, coefficients = coef)
})

results_df <- bind_rows(lapply(results, function(x){
  data.frame(
    formula = x$formula,
    AICc = x$AICc,
    R2 = x$R2
  )
}))

# Ordenar por AICc
results_df$deltaAICc = results_df$AICc - min(results_df$AICc)
results_df$weight = exp(-0.5 * results_df$deltaAICc)
results_df$weight = results_df$weight / sum(results_df$weight)
results_df <- results_df %>% arrange(deltaAICc)
results_df[1:10,]

write.csv(results_df, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_PGLS.csv")

best_pla <- pgls(normas ~ size, data = comp)
summary(best_pla)

plot(fitted(best_pla), residuals(best_pla),
     xlab = "Valores ajustados",
     ylab = "Resíduos",
     main = "Platyrrhini: resíduos vs ajustados")
abline(h = 0, col = "red", lty = 2)

hist(residuals(best_pla), main="Platyrrhini: Histograma de resíduos", xlab="Resíduos")
qqnorm(residuals(best_pla)); qqline(residuals(best_pla), col="red")

shapiro.test(residuals(best_pla))

res_pla <- residuals(best_pla)
lambda_pla <- phylosig(tree_p, res_pla, method="lambda")
lambda_pla

#################################################################################
#################################################################################
### EFEITO DA PARVORDER
df <- data.frame(
  species = mds_f$especies,
  parvorder = factor(mds_f$dados.PARVORDER),
  size = evolvas_f$Size,
  integration = evolvas_f$Integration,
  dimorfism = evolvas_f$Dimorfism,
  #cvn = nicho_f$CVn_mean,
  #clima = nicho_f$climacv_mean,
  social = mds_f$dados.SOCIAL_ORGANIZATION,
  mating = mds_f$dados.MATING_SYSTEM,
  prop = mds_f$dados.PROP_MALES_FEMALES,
  dominance = mds_f$dados.DOMINANCE,
  aggression = mds_f$dados.AGGRESSION,
  nmds1 = mds_f$scores.fit....1.,
  nmds2 = mds_f$scores.fit....2.,
  normas = align_f$normas
)

# Reordenar para bater com a árvore
df <- df[match(tree$tip.label, df$species), ]
rownames(df) <- df$species

# Preparar dados comparativos para PGLS
comp <- comparative.data(tree, df, names.col = species, vcv = TRUE)

model_basic <- pgls(normas ~ size + integration + nmds1 + nmds2 + parvorder, data = comp)
summary(model_basic)

model_inter <- pgls(normas ~ (size + integration + nmds1 + nmds2) * parvorder, data = comp)
summary(model_inter)

# Número de parâmetros
k1 <- attr(logLik(model_basic), "df")
k2 <- attr(logLik(model_inter), "df")
n <- nrow(df)

# AICc
AICc1 <- -2*logLik(model_basic) + 2*k1 + (2*k1*(k1+1))/(n - k1 - 1)
AICc2 <- -2*logLik(model_inter) + 2*k2 + (2*k2*(k2+1))/(n - k2 - 1)

delta <- c(AICc1, AICc2) - min(c(AICc1, AICc2))
weights <- exp(-0.5*delta) / sum(exp(-0.5*delta))

delta
weights

# Normalizar cvn para alpha se quiser
#df$cvn_scaled <- rescale(df$cvn, to = c(0.3, 1))

# Previsões do modelo final (pode escolher basic ou inter)
df$predicted <- as.numeric(predict(model_basic, newdata = df))

# Cores
mycols <- c("Catarrhini" = "#1b9e77", "Platyrrhini" = "#7570b3")

p1 <- ggplot(df, aes(x = size, y = normas)) +
  geom_text(aes(label = species, color = parvorder, size = integration, alpha = cvn)) +
  scale_size_continuous(name = "Integration") +
  scale_alpha_continuous(name = "CVn") +
  scale_color_manual(values = mycols) +
  theme_classic() +
  labs(x = "Average Skull Size", y = "Sexual Dimorphism", color = "Parvorder") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )

p1
