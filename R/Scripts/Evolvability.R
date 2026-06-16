# Evolvability of sexual dimorphism

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(ggplot2)
library(openxlsx)
library(evolqg)
library(evolvability)
library(ape)

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

Evolvability <- function(cov.matrix, beta.mat = NULL, iterations = 1000, size_skull) {
  # Funções auxiliares internas
  Norm <- function(x) {
    sqrt(sum(x * x))
  }
  
  Normalize <- function(x) {
    x / Norm(x)
  }
  
  # Número de traços
  num.traits <- dim(cov.matrix)[1]
  
  # Se beta.mat não for fornecida, sorteia vetores aleatórios e normaliza
  if (is.null(beta.mat)) {
    beta.mat <- array(rnorm(num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply(beta.mat, 2, Normalize)
  }
  
  # Calcula evolvabilidade (resposta)
  respostas <- diag(t(beta.mat) %*% cov.matrix %*% beta.mat)
  
  # Métricas
  respostas_pure <- mean(respostas)
  respostas_normal <- mean(respostas / size_skull)
  icv <- sd(respostas) / mean(respostas)
  icv_normal <- icv / size_skull
  
  # Retorna lista com resultados
  return(list(
    respostas = respostas,
    respostas_normal = respostas_normal,
    icv = icv,
    icv_normal = icv_normal,
    respostas_pure = respostas_pure
  ))
}

# read all VCV matrices
setwd("~/Dropbox/Doc/Data/vcv/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement = "", temp)
vcv$Lagothrix_lagothricha <- as.data.frame(lapply(vcv$Lagothrix_lagothricha, function(x) as.numeric(trimws(x))))
vcv$Cacajao_calvus = read.csv("~/Dropbox/Doc/Data/p_vcv_gabriel/Cacajao_calvus.csv", header = FALSE, sep = ";", dec = ",")
vcv$Cacajao_calvus <- as.data.frame(lapply(vcv$Cacajao_calvus, function(x) as.numeric(trimws(x))))

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts")
#
matings = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")
align = read.table("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv", header = TRUE, sep = ",")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/ancestrals_averages_PCS_autovalues_primates.RDS")

# Define species based on match between vcv and matings
commom <- intersect(names(vcv), matings$especies)

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
tree$tip.label[tree$tip.label == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"
species = matings$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))

# cut some vcv
vcv <- vcv[names(vcv) %in% matings$especies]
tree = drop.tip(tree, setdiff(matings$especies, names(vcv)))

all_cov_matrices = PhyloW(tree, vcv)
ancestral = getMRCA(tree, tree$tip.label)

results <- list()  # lista vazia
for(i in seq_along(matings$especies)){
  sp <- matings$especies[i]
  cat("Rodando:", sp, "\n")   # DEBUG: imprime no console
  
  medidas <- medias$ByTrait_Averages[[sp]]

  # checagens
  if (is.null(medidas)) { 
    cat("   Sem medidas para:", sp, "\n")
    next 
  }

  
  # se não tem medidas, pula
  if(is.null(medidas)){
    next
  }
  
  # se não está no vcv, pula
  if(!(sp %in% names(vcv))){
    next
  }
  
  # cálculos
  size <- (geomean(medidas$Machos) + geomean(medidas$Fêmeas)) / 2

  # calcula médias dos traits
  medias_por_trait <- (medidas$Machos + medidas$Fêmeas) / 2
  medias_por_trait <- as.numeric(medias_por_trait)
  
  # padronizar matriz
  covar <- vcv[[sp]]
  covar <- as.matrix(covar)                          # garante matriz numérica
  mean_prod <- outer(medias_por_trait, medias_por_trait, "*") 
  mat_scaled <- covar / mean_prod   # padroniza tudo
  diag(mat_scaled) <- diag(covar) / (medias_por_trait ^ 2)   
  covar <- mat_scaled
  
  if (is.null(covar)) { 
    cat("   Sem covar no vcv para:", sp, "\n")
    next 
  }
  
  rownames(covar) <- colnames(covar) <- NULL
  
  # normalizações
  dimor <- medidas$Machos - medidas$Fêmeas
  dimor_norm <- dimor / sqrt(sum(dimor ^ 2))
  
  # evolvability
  resposta <- as.matrix(covar) %*% as.vector(dimor_norm)
  evolv <- sum((resposta * dimor_norm) / (size ^ 2))
  
  evolv_media = unname(MeanMatrixStatistics(covar)[6])
  evolv_var = var(Evolvability(covar, size_skull = size)$respostas)
  
  # avg_evolvability OK
  avg_evolvability = sum(diag(covar)) / ncol(covar)
  
  # std evolvability OK
  std_evolv <- sqrt(evolv) / (size * size)
  
  ev <- eigen(covar, symmetric = TRUE)
  ratio_peak_mean <- ev$values[1] / mean(ev$values)
  
  proj <- t(ev$vectors) %*% dimor_norm
  contrib <- (proj^2) / sum(proj^2)
  
  maior_lambda = eigen(covar)$values[1]
  
  # integration OK
  rownames(covar) <- colnames(covar) <- NULL
  integ <- CalcEigenVar(covar)
  integ2 <- CalcR2(cov2cor(covar))
  integ_media = 1 - unname(MeanMatrixStatistics(covar))[8] # mesmos valores para cor e cov

  # salva numa linha
  results[[sp]] <- data.frame(
    Species = sp,
    Evolvability = evolv,
    Standart_Evolvability = std_evolv,
    Dimorfism = norma(dimor),
    Size = size,
    Integration = integ,
    Average_Evolvability = avg_evolvability,
    Ratio_Peak_Mean = ratio_peak_mean,
    Contribution = contrib[1],
    Maior_Lambda = maior_lambda,
    Evolvability_Media = evolv_media,
    Evolvability_Variance = evolv_var,
    Integration_Media = integ_media,
    Integration_R2 = integ2
  )
}

# junta tudo num único data.frame
evolvas <- do.call(rbind, results)
rownames(evolvas) <- NULL

# save results
#saveRDS(evolvas, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
#evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")

evolvas <- evolvas[match(tree$tip.label, evolvas$Species), ]
evolvas$genus <- factor(evolvas$Species, levels = tree$tip.label)

# dimorfismo
plot(evolvas$Integration, evolvas$Evolvability)

# 1 - autonomy media
plot(evolvas$Integration_Media, evolvas$Evolvability_Media)
plot(evolvas$Integration_Media, evolvas$Evolvability_Variance)

#EigenVar
plot(evolvas$Integration, evolvas$Evolvability_Variance)
plot(evolvas$Integration, evolvas$Evolvability_Media)

p <- ggplot(
  evolvas,
  aes(
    x = Dimorfism,
    y = Evolvability
  )
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
    plot.title = element_text(
      size = 28,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 20)
    ),
    axis.title = element_text(
      size = 28,
      face = "bold"
    ),
    axis.text = element_text(
      size = 24,
      face = "bold",
      colour = "black"
    ),
    legend.title = element_text(
      size = 24,
      face = "bold"
    ),
    legend.text = element_text(
      size = 24,
      face = "bold"
    )
  ) +
  
  labs(
    x = "Global Integration",
    y = "Average Evolvability"
  )

p

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Urgente6_Evolvability_Integration_Pmax_Error_NonError.png",
  plot = p1,
  width = 14,
  height = 7,
  dpi = 300
)

cor_matrix <- cor(
  evolvas[, c("Integration", "Integration_Media", "Integration_R2")],
  use = "complete.obs",
  method = "pearson"
)

cor_matrix

# CORRELACAO ENTRE INTEGRACAO E 1 - AUTONOMIA  OK
# CONTROLAR A MATRIZ DE COVARIANCIA


p1 = ggplot(evolvas, aes(x = Species, y = Evolvability)) +
     geom_col() +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_discrete(drop = TRUE) +
     xlab(" Species") +
     ylab("Evolvability SD") +
     coord_flip() +
     theme_classic()
p1

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability_Species.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


cat("i =", i, 
    "| species[i] =", species[i], 
    "| vcv name =", names(vcv)[i], "\n")

#summary(evolvas)

evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Dimorfism, y = Evolvability), size = 2) +
  #geom_text(aes(x = Dimorfism, y = Evolvability, label = species)) +
  #geom_smooth(aes(x = Dimorfism, y = Evolvability), method = "lm") +
  xlab("Magnitude of the Sexual Dimorfism") +
  ylab("Evolvability on the Direction of Sexual Dimorphism") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14)
  )

evolv_plot

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability_Dimorphism.png", plot = evolv_plot,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução
