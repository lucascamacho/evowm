# PGLS POR GENERO
# CATARRHINI, PLATYRRHINI, HAPLORRHINI
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

library(slouch)
library(ape)
library(dplyr)
library(caper)
library(phytools)
library(ggplot2)
library(ggrepel)
library(scales)  # para escala de transparência
library(gt)
library(car)

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")
#nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
  #nicho$genus,
  align$matings.genus,
  evolvas$genus,
  mds$genus
))

#nicho_f   <- nicho[nicho$genus %in% sp_comuns, ]
align_f   <- align[align$matings.genus %in% sp_comuns, ]
evolvas_f <- evolvas[evolvas$genus %in% sp_comuns, ]
mds_f     <- mds[mds$genus %in% sp_comuns, ]

# read all vcv matrices
setwd("~/Dropbox/Doc/Data/genus_vcv/")
temp = list.files(pattern = "*.txt")
vcv = lapply(temp, read.csv, header = FALSE, dec = ".", sep = "\t")
names(vcv)  = gsub(".csv", replacement = "", temp)

wrong <- which(sapply(vcv, ncol) == 1)

for(i in wrong){
  
  file <- temp[i]
  
  # tenta outras combinações
  mat <- read.table(file, header = FALSE, dec = ".")
  
  vcv[[i]] <- as.matrix(mat)
}
vcv <- lapply(vcv, as.matrix)
names(vcv) <- gsub("\\.txt$", "", temp)

# read and plot phylo tree
commom <- intersect(names(vcv), mds_f$genus)
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
genus <- sub("_.*", "", tree$tip.label)
keep <- genus %in% commom
tree2 <- drop.tip(tree, tree$tip.label[!keep])
tree <- tree2
# extrair generos
genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
unique_genera <- unique(genera)
representatives <- c()
for(g in unique_genera){
  
  spp <- tree$tip.label[genera == g]
  
  if(length(spp) == 1){
    
    representatives[g] <- spp
    
  } else {
    
    # pegar comprimento do ramo terminal
    edges <- match(spp, tree$tip.label)
    terminal_lengths <- tree$edge.length[match(edges, tree$edge[,2])]
    
    # escolher o menor (mais recente)
    representatives[g] <- spp[which.max(terminal_lengths)]
  }
}

# remover todas as outras espécies
tree_genus <- drop.tip(tree, setdiff(tree$tip.label, representatives))

# renomear para os gêneros
tree_genus$tip.label <- names(representatives)

#
# media_genero <- nicho %>%
#   dplyr::group_by(.data$genus) %>%
#   dplyr::summarise(
#   climacv_mean = mean(.data$climacv, na.rm = TRUE),
#   CVn_mean = mean(.data$CVn, na.rm = TRUE)
#   )

#nicho_f <- media_genero[match(tree_genus$tip.label, media_genero$genus), ]
evolvas_f <- evolvas_f[match(tree_genus$tip.label, evolvas_f$genus), ]
align_f <- align_f[match(tree_genus$tip.label, align_f$matings.genus), ]
mds_f <- mds_f[match(tree_genus$tip.label, mds_f$genus), ]

################################################################################
# A PARTIR DAQUI
# Catarrhini #1:21 sem clima e cvn, 1:16 com clima e cvn
align_c = align_f[1:21,]
evolvas_c = evolvas_f[1:21,]
mds_c = mds_f[1:21,]
tree_c = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_c$matings.genus))
#nicho_c = nicho_f[1:16,]

#
df <- data.frame(
  genus = mds_c$genus,
  parvorder = mds_c$dados.PARVORDER,
  social = mds_c$dados.SOCIAL_ORGANIZATION,
  mating = mds_c$dados.MATING_SYSTEM,
  prop = mds_c$dados.PROP_MALES_FEMALES,
  dominance = mds_c$dados.DOMINANCE,
  aggression = mds_c$dados.AGGRESSION,
  all_male = mds_c$dados.ALL_MALE_GROUPS,
  furtive = mds_c$dados.FURTIVE_COPULATION,
  dimorfism = evolvas_c$Dimorfism,
  integration = evolvas_c$Integration,
  size = evolvas_c$Size,
  nmds1 = mds_c$vegan..scores.fit..sites...1.,
  nmds2 = mds_c$vegan..scores.fit..sites...2.,
  normas = align_c$normas
  #clima = nicho_c$climacv_mean,
  #cvn = nicho_c$CVn_mean
)

df <- df[match(tree_c$tip.label, df$genus), ]

#predictors <- c("size", "nmds1", "nmds2", "clima", "cvn", "mating", "social", "prop", "dominance", "integration")
predictors <- c("size", "nmds1", "nmds2", "integration")
response <- "normas"

rownames(df) = df$genus  # ou use df$genus se sua árvore for por gênero

# Cria o objeto comparative.data
comp <- comparative.data(
  phy = tree_c,       # sua árvore
  data = df,          # seu dataframe
  names.col = "genus", # coluna com os nomes que batem com a árvore
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
  
  # lambda dos resíduos
  res <- residuals(fit)
  lambda <- phylosig(tree_c, res, method = "lambda")$lambda
  
  list(
    formula = f,
    fit = fit,
    AICc = aicc,
    R2 = r2,
    lambda = lambda,
    coefficients = coef
  )
})

results_df <- bind_rows(lapply(results, function(x){
  data.frame(
    formula = x$formula,
    AICc = x$AICc,
    R2 = x$R2,
    lambda = x$lambda
  )
}))

# Ordenar por AICc
results_df$deltaAICc = results_df$AICc - min(results_df$AICc)
results_df$weight = exp(-0.5 * results_df$deltaAICc)
results_df$weight = results_df$weight / sum(results_df$weight)
results_df <- results_df %>% arrange(deltaAICc)
results_df[1:10,]

#write.csv(results_df, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/catarrhini_PGLS.csv")

#best_cat <- pgls(integration ~ nmds1 + clima + cvn, data = comp)
best_cat <- pgls(normas ~ size + integration, data = comp)


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
# Platyrrhini # 22:37 sem clima e cvn, 17:28 com clima e cvn
align_p = align_f[22:37,]
evolvas_p = evolvas_f[22:37,]
mds_p = mds_f[22:37,]
tree_p = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_p$matings.genus))
#nicho_p = nicho_f[17:28,]

#
df <- data.frame(
  genus = mds_p$genus,
  parvorder = mds_p$dados.PARVORDER,
  social = mds_p$dados.SOCIAL_ORGANIZATION,
  mating = mds_p$dados.MATING_SYSTEM,
  prop = mds_p$dados.PROP_MALES_FEMALES,
  dominance = mds_p$dados.DOMINANCE,
  aggression = mds_p$dados.AGGRESSION,
  furtive = mds_p$dados.FURTIVE_COPULATION,
  infanticide = mds_p$dados.INFANTICIDE,
  dimorfism = evolvas_p$Dimorfism,
  integration = evolvas_p$Integration,
  size = evolvas_p$Size,
  nmds1 = mds_p$vegan..scores.fit..sites...1.,
  nmds2 = mds_p$vegan..scores.fit..sites...2.,
  normas = align_p$normas
  #clima = nicho_p$climacv_mean,
  #cvn = nicho_p$CVn_mean
)

df <- df[match(tree_p$tip.label, df$genus), ]

#predictors <- c("size", "nmds1", "nmds2",  "integration", "furtive", "prop", "social", "mating", "clima", "cvn")
predictors <- c("size", "nmds1", "nmds2", "integration")
response <- "normas"

rownames(df) = df$genus  # ou use df$genus se sua árvore for por gênero

# Cria o objeto comparative.data
comp <- comparative.data(
  phy = tree_p,       # sua árvore
  data = df,          # seu dataframe
  names.col = "genus", # coluna com os nomes que batem com a árvore
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
  
  # lambda dos resíduos
  res <- residuals(fit)
  lambda <- phylosig(tree_p, res, method = "lambda")$lambda
  
  list(
    formula = f,
    fit = fit,
    AICc = aicc,
    R2 = r2,
    lambda = lambda,
    coefficients = coef
  )
})

results_df <- bind_rows(lapply(results, function(x){
  data.frame(
    formula = x$formula,
    AICc = x$AICc,
    R2 = x$R2,
    lambda = x$lambda
  )
}))

# Ordenar por AICc
results_df$deltaAICc = results_df$AICc - min(results_df$AICc)
results_df$weight = exp(-0.5 * results_df$deltaAICc)
results_df$weight = results_df$weight / sum(results_df$weight)
results_df <- results_df %>% arrange(deltaAICc)
results_df[1:10,]

#write.csv(results_df, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/platyrrhini_PGLS.csv")

#best_pla <- pgls(integration ~ normas + nmds1 + prop + dominance, data = comp)
best_pla <- pgls(normas ~ size + nmds2 + integration, data = comp)

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
# PLOT TABLE 
catarrhini_df <- data.frame(
  Clade = "Catarrhini",
  Model = c(
    "Size + Integration",
    "Size + NMDS1 + Integration",
    "Size + NMDS2 + Integration"
  ),
  AICc = c(14.10283, 15.57243, 16.71276),
  R2 = c(0.7324028, 0.7376805, 0.7230423),
  lambda = c(0.9538401, 0.9999267, 0.9285472),
  deltaAICc = c(0.000000, 1.469599, 2.609935),
  weight = c(0.4500330, 0.2158366, 0.1220405)
)

platyrrhini_df <- data.frame(
  Clade = "Platyrrhini",
  Model = "Size + NMDS2 + Integration",
  AICc = -8.529196,
  R2 = 0.7206486,
  lambda = 1.311091,
  deltaAICc = 0.000000,
  weight = 0.5467892029
)

catarrhini_df <- catarrhini_df %>%
  filter(deltaAICc < 2)

platyrrhini_df <- platyrrhini_df %>%
  filter(deltaAICc < 2)

results_pgls <- bind_rows(catarrhini_df, platyrrhini_df)

results_pgls <- results_pgls %>%
  mutate(
    Clade = case_when(
      Clade == "Catarrhini" ~ "Catarrhini (N = 27)",
      Clade == "Platyrrhini" ~ "Platyrrhini (N = 16)"
    )
  )

gt_pgls <- results_pgls %>%
  
  gt(groupname_col = "Clade") %>%
  
  cols_label(
    Model = "Model",
    AICc = "AICc",
    R2 = "R²",
    lambda = "λ",
    deltaAICc = "ΔAICc",
    weight = "Weight"
  ) %>%
  
  fmt_number(
    columns = c(AICc, R2, lambda, deltaAICc, weight),
    decimals = 3
  ) %>%
  
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  
  tab_header(
    title = md("**PGLS model selection (ΔAICc < 2)**")
  ) %>%
  
  tab_source_note(
    source_note = md(
      "Only models with ΔAICc < 2 are shown. AICc weights represent relative support among candidate models."
    )
  )


gt_pgls

gtsave(
  gt_pgls,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/pgls_table.png",
  vwidth = 1600,
  vheight = 1000,
  zoom = 2
)

#################################################################################
### EFEITO DA PARVORDER
df <- data.frame(
  genus = mds_f$genus,
  parvorder = factor(mds_f$dados.PARVORDER),
  size = evolvas_f$Size,
  integration = evolvas_f$Integration,
  dimorfism = evolvas_f$Dimorfism,
  nmds1 = mds_f$vegan..scores.fit..sites...1.,
  nmds2 = mds_f$vegan..scores.fit..sites...2.,
  social = mds_f$dados.SOCIAL_ORGANIZATION,
  mating = mds_f$dados.MATING_SYSTEM,
  prop = mds_f$dados.PROP_MALES_FEMALES,
  dominance = mds_f$dados.DOMINANCE,
  aggression = mds_f$dados.AGGRESSION,
  normas = align_f$normas
)

# Reordenar para bater com a árvore
df <- df[match(tree_genus$tip.label, df$genus), ]
rownames(df) <- df$genus

vars <- df[, c("size", "integration", "nmds1", "nmds2")]

cor_matrix <- cor(vars, method = "pearson", use = "complete.obs")
round(cor_matrix, 2)


# Preparar dados comparativos para PGLS
comp <- comparative.data(tree_genus, df, names.col = genus, vcv = TRUE)

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

res_haplo <- residuals(model_basic)
lambda_haplo <- phylosig(tree_genus, res_haplo, method="lambda")
lambda_haplo


# Previsões do modelo final (pode escolher basic ou inter)
df$predicted <- as.numeric(predict(model_basic, newdata = df))

df <- df[order(df$integration), ]

r2 <- summary(model_basic)$r.squared

fstat <- summary(model_basic)$fstatistic
p_model <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

label <- paste0("R² = ", round(r2, 3),
                "\nP = ", format.pval(p_model, digits = 3, eps = 1e-3))
# Cores
mycols <- c("Catarrhini" = "#1b9e77", "Platyrrhini" = "#7570b3")

p1 <- ggplot(df, aes(x = integration, y = normas)) +
  geom_smooth(aes(y = predicted),
              method = "lm",
              se = TRUE,
              color = "black") +
  #geom_point(aes(color = parvorder, size = size, alpha = nmds1)) +
  geom_point(aes(color = parvorder, size = size)) +
  annotate("text",
           x = Inf, y = Inf,
           label = label,
           hjust = 1.1, vjust = 1.5,
           size = 5) +
  scale_size_continuous(name = "Size") +
  #scale_alpha_continuous(range = c(0.1, 1), limits = c(-0.1, 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6, col = "red") +
  scale_color_manual(values = mycols) +
  theme_classic() +
  labs(x = "Integration Index", y = "Sexual Dimorphism", color = "Parvorder") +
  #labs(x = "Integration Index", y = "Sexual Dimorphism", color = "Parvorder", alpha = "NMDS1") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

p1

# Salva o gráfico em alta resolução
ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/BestModel_Haplorrhini.png",
  plot = p1,
  width = 12,
  height = 7,
  dpi = 300
)

#
coefs <- summary(model_basic)$coefficients

tab <- data.frame(
  term = rownames(coefs),
  estimate = coefs[,1],
  std.error = coefs[,2],
  statistic = coefs[,3],
  p.value = coefs[,4]
)

#-----------------------------
# 2. Renomear preditores
#-----------------------------
tab$term <- dplyr::recode(tab$term,
                          "(Intercept)" = "Intercept",
                          "size" = "Size",
                          "integration" = "Integration",
                          "nmds1" = "NMDS1",
                          "nmds2" = "NMDS2",
                          "parvorderPlatyrrhini" = "Parvorder"
)

#-----------------------------
# 3. Formatação (sem sig styling)
#-----------------------------
tab <- tab %>%
  mutate(
    est_fmt = round(estimate, 3),
    se_fmt = round(std.error, 3),
    p_fmt = ifelse(
      p.value < 0.001,
      "<0.001",
      as.character(round(p.value, 3))
    )
  )

#-----------------------------
# 4. Tabela gt
#-----------------------------
gt_parvo = gt(tab %>% dplyr::select(term, est_fmt, se_fmt, statistic, p_fmt)) %>%
  
  cols_label(
    term = "Predictor",
    est_fmt = "Estimate",
    se_fmt = "SE",
    statistic = "t",
    p_fmt = "p"
  ) %>%
  
  fmt_number(
    columns = c(est_fmt, se_fmt, statistic),
    decimals = 3
  ) %>%
  
  tab_header(
    title = md("**PGLS model including parvorder effect**")
  ) %>%
  
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(everything())
  ) %>%
  
  tab_options(
    table.align = "center"
  )

gtsave(
  gt_parvo,
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/pgls_parvorder.png",
  vwidth = 1600,
  vheight = 500,
  zoom = 2
)
