# estimating k1 and k2
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability")

library(ape)
library(ggplot2)
library(evolvability)
library(evolqg)
library(tidyverse)

# 
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
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
species <- species[1:49]

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
  integ <- CalcEigenVar(as.matrix(covar))
 
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

saveRDS(error, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Error.RDS")
#readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Error.RDS")

error_long <- error %>%
  dplyr::select(
    species,
    Dimorfism,
    Evolvability,
    Evolvability_error_mean
  ) %>%
  pivot_longer(
    cols = c(Evolvability, Evolvability_error_mean),
    names_to = "Type",
    values_to = "Evolvability_value"
  ) %>%
  mutate(
    Type = recode(
      Type,
      "Evolvability" = "Observed",
      "Evolvability_error_mean" = "Error-corrected"
    )
  )

ggplot(error_long,
       aes(x = Evolvability_value,
           y = Dimorfism,
           color = Type)) +
  
  # linhas ligando os pares
  geom_segment(
    data = error,
    aes(
      x = Evolvability,
      xend = Evolvability_error_mean,
      y = Dimorfism,
      yend = Dimorfism
    ),
    inherit.aes = FALSE,
    color = "grey75",
    linewidth = 0.8,
    alpha = 0.7
  ) +
  
  # pontos com jitter
  geom_jitter(
    width = 0.0006,
    height = 0.03,
    size = 3.5,
    alpha = 0.9
  ) +
  
  scale_color_manual(values = c(
    "Observed" = "#1b9e77",
    "Error-corrected" = "#d95f02"
  )) +
  
  labs(
    x = "Evolvability",
    y = "Sexual dimorphism",
    color = NULL,
    title = "Observed vs error-corrected evolvability"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

