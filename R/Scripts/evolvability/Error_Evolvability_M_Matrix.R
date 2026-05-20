# MICROSCRIBE AND POLHEMUS MEDE A MESMA COISA COM DIFENRENTES SIZES
# M MATRIX E GRANDE 
# LINEAR MODEL FOR SIZE AND PESSOA QUE MEDIU USANDO TODOS OS VALORES
# COEF OF DETERMINATION (1 - r2) OK
# MOSTRAR O PLOT ORIGINAL COM SIZE NO BANNER OK

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
# 
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

#
autonomia <- function(R) {
  diag_s <- numeric(ncol(R))
  
  for (i in 1:ncol(R)) {
    r2 <- summary(lm(R[, i] ~ R[, -i]))$r.squared
    diag_s[i] <- r2
  }
  
  mean(1 - diag_s)
}

integration_index <- function(R) {
  a <- autonomia(R)
  1 - a
}

#
prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

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
  
  #
  dimor <- medidas$Machos - medidas$Fêmeas
  
  #
  covar <- vcv[[sp]]
  covar <- as.matrix(covar)
  
  if (is.null(covar)) { 
    cat("   Sem covar no vcv para:", sp, "\n")
    next 
  }
  
  # 
  eig <- eigen(covar)
  D <- diag(eig$values)
  V <- eig$vectors
  D2 <- D
  D2[1,1] <- 0
  covar <- V %*% D2 %*% t(V)
  
  #
  dimor_norm <- dimor / sqrt(sum(dimor ^ 2))
  
  # 
  evolv <- as.numeric(t(dimor_norm) %*% covar %*% dimor_norm)
  
  # 
  rownames(covar) <- colnames(covar) <- NULL
  correl <- cov2cor(covar)
  
  integ <- CalcEigenVar(as.matrix(correl))
  #integ <- unname(1 - evolvability::evolvabilityBeta(correl, dimor_norm)$a)
  #integ <- integration_index(correl)
  
  #
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
    
    B_err <- B + e
    
    evolv_errors[k] <- as.numeric(
      t(B_err) %*% P %*% B_err /
        (t(B_err) %*% B_err)
    )
  }
  
  #
  results[[sp]] <- data.frame(
    species = sp,
    
    Evolvability = evolv,
    
    Evolvability_error_mean = mean(evolv_errors),
    Evolvability_error_sd = sd(evolv_errors),
    
    Evolvability_error_q025 = quantile(evolv_errors, 0.025),
    Evolvability_error_q975 = quantile(evolv_errors, 0.975),
    
    Dimorfism = norma(dimor),
    Size = size,
    Integration = integ
  )
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
    Integration
  ) %>%
  rename(
    Evolvability_value = Evolvability
  )

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
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24, face = "bold")
  ) +
  
  labs(
    x = "Sexual Dimorphism",
    y = "Evolvability"
  )

p1

# # Salva o gráfico em alta resolução
   ggsave(
     "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Dimorphism_Error_NonError.png",
      plot = p1,
      width = 12,
      height = 7,
      dpi = 300
   )

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
      x = "Integration Index",
      y = "Evolvability"
    )
  
  p2

  ggsave(
    "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Integration_Error_NonError.png",
    plot = p2,
    width = 12,
    height = 7,
    dpi = 300
  )
  
