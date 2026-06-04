# estimating k1 and k2
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability")

library(ape)
library(ggplot2)
library(evolvability)
library(evolqg)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggrepel)
library(ggpmisc)
library(MASS)

# 
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

#
prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))

#
setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
temp <- list.files(pattern = "*.csv")
vcv <- lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
names(vcv) <- gsub(".csv", replacement= "", temp)

#
error_samples <- readRDS("~/Dropbox/Doc/Data/error_samples.RDS")
names(error_samples) <- gsub(" ", "_", names(error_samples))

#
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")

#
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = read.nexus(filename)
species = names(vcv)
species = append(species, "Homo_sapiens")
tree = drop.tip(tree, setdiff(tree$tip.label, species))

species <- species[match(tree$tip.label, species)]

#
results <- list()  # lista vazia

for(i in seq_along(species)){
  sp <- species[i]
  cat("Rodando:", sp, "\n")   #
  
  medidas <- medias$ByTrait_Averages[[sp]]
  
  # 
  if (is.null(medidas)) { 
    cat("   Sem medidas para:", sp, "\n")
    next 
  }

  # 
  if(is.null(medidas)){
    next
  }
  
  # 
  if(!(sp %in% names(vcv))){
    next
  }
  
  #
  if(!(sp %in% names(error_samples))){
    cat("   Sem error_samples para:", sp, "\n")
    next
  }
  
  # 
  size <- (geomean(medidas$Machos) + geomean(medidas$Fêmeas)) / 2
  size_vec <- (medidas$Machos + medidas$Fêmeas) / 2
  
  
  #
  dimor <- medidas$Machos - medidas$Fêmeas
  
  #
  covar <- vcv[[sp]]
  covar <- as.matrix(covar)
  
  if (is.null(covar)) { 
    cat("   Sem covar no vcv para:", sp, "\n")
    next 
  }
  
   eig <- eigen(covar)
   D <- diag(eig$values)
   V <- eig$vectors
   D2 <- D
   D2[1,1] <- min(eig$values[-1]) * 1e-8
   covar <- V %*% D2 %*% t(V)
  
  #
  dimor_norm <- dimor / sqrt(sum(dimor ^ 2))
  
  # 
  evolv <- as.numeric(t(dimor_norm) %*% covar %*% dimor_norm)
  
  eig <- eigen(covar)
  pmax <- eig$vectors[, 1]
  pmax <- pmax / sqrt(sum(pmax^2))  # normaliza
  e_pmax <- as.numeric(t(pmax) %*% covar %*% pmax)
  
  # 
  rownames(covar) <- colnames(covar) <- NULL
  correl <- cov2cor(covar)
  
  #integ <- CalcEigenVar(as.matrix(correl))
  integ <- unname(1 - evolvability::evolvabilityBeta(correl, dimor_norm)$a)
  #integ <- unname(evolvability::evolvabilityBeta(correl, dimor_norm)$i)
  
  
  mat_errors <- error_samples[[sp]]
  
  measure_cols <- c(
    "ISPM","ISNSL","ISPNS","PMZS","PMZI","PMMT",
    "NSLNA","NSLZS","NSLZI","NABR","NAFM","NAPNS",
    "BRPT","BRAPET","PTFM","PTAPET","PTBA","PTEAM",
    "PTZYGO","FMZS","FMMT","ZSZI","ZIMT",
    "ZIZYGO","ZITSP","MTPNS","PNSAPET","APETBA",
    "APETTS","BAEAM","EAMZYGO","ZYGOTSP","LDAS",
    "BRLD","OPILD","PTAS","JPAS","BAOPI"
  )
  evolv_errors <- numeric(nrow(mat_errors))
  
  for(k in 1:nrow(mat_errors)){
    P <- covar
    B <- dimor_norm
    
    e <- as.numeric(mat_errors[k, measure_cols])
    
    B_err_raw <- B + e
    
    # 2. Re-unitiza o vetor com erro (FORÇA o comprimento a ser 1 novamente)
    B_err <- B_err_raw / sqrt(sum(B_err_raw ^ 2))
    
    # 3. Calcula a evolvabilidade pura (sem precisar dividir por nada embaixo)
    evolv_errors[k] <- as.numeric(t(B_err) %*% P %*% B_err)
    
  }
  
  p1 <- eig$vectors[,1]
  p1 <- p1 / sqrt(sum(p1^2))
  
  a_pmax <- as.numeric(t(dimor_norm) %*% p1) * p1
  
  e_dimor_pmax <- as.numeric(t(a_pmax) %*% covar %*% a_pmax)
  
  
  size_vec_norm <- size_vec/ sqrt(sum(size_vec^2))
  
  e_size <- as.numeric(t(size_vec_norm) %*% covar %*% size_vec_norm)
  
  #options(warn = -1)
  results[[sp]] <- data.frame(
    species = sp,
    
    Evolvability = evolv,
    Mean_Evolvability = sum(diag(covar)) / ncol(covar),
    
    Conditional_Evolvability = 1 / (t(dimor_norm) %*% solve(covar) %*% dimor_norm),
    Mean_Conditional_Evolvability = (length(eigen(covar)$values) - 1) / sum(1 / eigen(covar)$values[-1]),
    
    Pmax = e_pmax,
    Evolvability_error_mean = mean(evolv_errors),
    Evolvability_error_sd = sd(evolv_errors),
    
    Evolvability_error_q025 = quantile(evolv_errors, 0.025),
    Evolvability_error_q975 = quantile(evolv_errors, 0.975),
    
    Dimorfism = norma(dimor),
    Size = size,
    
    Integration = integ,
    
    Evolv_Dimor_Pmax = e_dimor_pmax,
    Evolv_Size = e_size
  )
  #options(warn = 0)
} 

# junta tudo num único data.frame
error <- do.call(rbind, results)
rownames(error) <- NULL

#saveRDS(error, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Error_M_Matrix.RDS")
#error <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Error_M_Matrix.RDS")

error_long <- error %>%
  dplyr::select(
    species,
    Dimorfism,
    Evolvability,
    Evolvability_error_q025,
    Evolvability_error_q975,
    Mean_Evolvability,
    Conditional_Evolvability,
    Mean_Conditional_Evolvability,
    Pmax,
    Integration,
    Evolv_Dimor_Pmax,
    Evolv_Size
  ) %>%
  rename(
    Evolvability_value = Evolvability
  )

plot(error$Evolvability, error$Evolv_Dimor_Pmax)
summary(lm(error$Evolvability ~ error$Evolv_Dimor_Pmax))

p1 <- ggplot(
  error_long,
  aes(
    x = Dimorfism,
    y = Evolvability_value
  )
) +
  
  # error bars
  geom_errorbar(
    aes(
      ymin = Evolvability_error_q025,
      ymax = Evolvability_error_q975
    ),
    width = 0,
    alpha = 0.5
  ) +
  
  # pontos
  geom_quasirandom(
    width = 0.02,
    size = 6,
    alpha = 0.8
  ) +
  
  # regressão
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  # R²
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.x = "right",
    label.y = "top",
    size = 10
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24, face = "bold")
  ) +
  
  labs(
    title = "Positive Relation Between Sexual Dimorphism and Evolvability",
    x = "Sexual Dimorphism",
    y = "Evolvability"
  )

p1

# # Salva o gráfico em alta resolução
   # ggsave(
   #   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Dimorphism_Error_NonError.png",
   #    plot = p1,
   #    width = 16,
   #    height = 7,
   #    dpi = 300
   # )

p2 <- ggplot(
    error_long,
    aes(
      x = Integration,
      y = Evolvability_value
    )
  ) +
    
    geom_errorbar(
      aes(
        ymin = Evolvability_error_q025,
        ymax = Evolvability_error_q975
      ),
      width = 0,
      alpha = 0.5
    ) +
    
    geom_point(
      size = 6,
      alpha = 0.8
    ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE
    ) +
    
    stat_poly_eq(
      aes(label = paste(..rr.label..)),
      formula = y ~ x,
      parse = TRUE,
      label.x = "left",
      label.y = "top",
      size = 10
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 28, face = "bold"),
      axis.text = element_text(
        size = 24,
        face = "bold",
        colour = "black"
      ),
      legend.title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 24, face = "bold")
    ) +
    
    labs(
      title = "No Relation between Constraint and Evolvability",
      x = "Constraint Index",
      y = "Evolvability"
  )
  
p2

# ggsave(
#    "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Integration_Pmax_Error_NonError.png",
#    plot = p2,
#    width = 14,
#    height = 7,
#    dpi = 300
# )

plot_e <- bind_rows(
  
  error_long %>%
    transmute(
      species,
      Dimorfism,
      Value = Evolvability_value,
      Type = "Sexual Dimorphism Evolvability"
    ),
  
  error_long %>%
    transmute(
      species,
      Dimorfism,
      Value = Mean_Evolvability,
      Type = "Mean Evolvability"
    )
)  
  
p3 <- ggplot() +
  
  geom_segment(
    data = error_long,
    aes(
      x = Dimorfism,
      xend = Dimorfism,
      y = Mean_Evolvability,
      yend = Evolvability_value
    ),
    alpha = 0.1
  ) +
  
  geom_point(
    data = plot_e,
    aes(
      x = Dimorfism,
      y = Value,
      shape = Type,
      color = Type
    ),
    size = 6,
    alpha = 0.9
  ) +
  labs(
    x = "Sexual Dimorphism",
    y = "Evolvability"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(
      size = 24,
      face = "bold",
      colour = "black"
    ),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24, face = "bold"),
    legend.position = "none"
  )  

p3  

#ggsave(
#   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Mean_Dimorphism.png",
#   plot = p3,
#   width = 16,
#   height = 7,
#   dpi = 300
#)

plot_ce <- bind_rows(
  
  error_long %>%
    transmute(
      species,
      Dimorfism,
      Value = Conditional_Evolvability,
      Type = "Sexual Dimorphism Conditional Evolvability"
    ),
  
  error_long %>%
    transmute(
      species,
      Dimorfism,
      Value = Mean_Conditional_Evolvability,
      Type = "Mean conditional evolvability"
    )
)


p4 <- ggplot() +
  
  geom_segment(
    data = error_long,
    aes(
      x = Dimorfism,
      xend = Dimorfism,
      y = Mean_Conditional_Evolvability,
      yend = Conditional_Evolvability
    ),
    alpha = 0.1
  ) +
  
  geom_point(
    data = plot_ce,
    aes(
      x = Dimorfism,
      y = Value,
      shape = Type,
      color = Type
    ),
    size = 6,
    alpha = 0.9
  ) +
  labs(
    x = "Sexual Dimorphism",
    y = "Conditional Evolvability"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(
      size = 24,
      face = "bold",
      colour = "black"
    ),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24, face = "bold")
    #legend.position = "none"
  ) 

p4

#ggsave(
#  "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Conditional_Evolvability_Mean_Dimorphism.png",
#  plot = p4,
#  width = 16,
#  height = 7,
#  dpi = 300
#)

p5 = ggplot(error_long, aes(x = Dimorfism, y = Evolvability_value)) +
  
  geom_point(color = "#00BFC4", shape = 17, size = 6, alpha = 0.8) +
  
  #scale_y_log10() +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(
      size = 24,
      face = "bold",
      colour = "black"
    )
  ) +
  
  labs(
    x = "Magnitude of sexual dimorphism",
    y = "Evolvability"
  )

# Para visualizar:
p5

# ggsave(
#   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_semPmax_Dimorphism.png",
#   plot = p5,
#   width = 16,
#   height = 7,
#   dpi = 300
# )

p6 = ggplot(error_long, aes(x = Dimorfism, y = Conditional_Evolvability)) +
  
  geom_point(size = 6, alpha = 0.9) +
  
  #scale_y_log10() +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(
      size = 24,
      face = "bold",
      colour = "black"
    )
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  stat_poly_eq(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left",
    label.y = "top",
    size = 10
  ) +
  
  labs(
    x = "Magnitude of sexual dimorphism",
    y = "Conditional Evolvability (log scale)"
  )

# Para visualizar:
p6

# ggsave(
#   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Conditional_Evolvability_Dimorphism.png",
#   plot = p6,
#   width = 16,
#   height = 7,
#   dpi = 300
# )

