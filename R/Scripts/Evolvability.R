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
  
  if (is.null(covar)) { 
    cat("   Sem covar no vcv para:", sp, "\n")
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
  # tam_cra <- geomean(medidas$Machos)
  size <- (geomean(medidas$Machos) + geomean(medidas$Fêmeas)) / 2
  #size <- (geomean(log(medidas$Machos)) + geomean(log(medidas$Fêmeas))) / 2
  # size <- ((log(medidas$Machos[2]) + log(medidas$Fêmeas[2])) /2) +
  # ((log(medidas$Machos[7]) + log(medidas$Fêmeas[7])) /2) +
  # ((log(medidas$Machos[10]) + log(medidas$Fêmeas[10])) /2) +
  # ((log(medidas$Machos[34]) + log(medidas$Fêmeas[34])) /2)
  
  # dimor <- (medidas$Machos - medidas$Fêmeas) / size
  dimor <- medidas$Machos - medidas$Fêmeas
  # dimor <- log(medidas$Machos) - log(medidas$Fêmeas)
  # calcula médias dos traits
  #medias_por_trait <- (medidas$Machos + medidas$Fêmeas) / 2
  #medias_por_trait <- as.numeric(medias_por_trait)   # garante vetor numérico
  
  # padronizar matriz
  covar <- vcv[[sp]]
  covar <- as.matrix(covar)                          # garante matriz numérica
  #mean_prod <- outer(medias_por_trait, medias_por_trait, "*")   # matriz com produtos das médias
  #mat_scaled <- covar / mean_prod   # padroniza tudo
  #diag(mat_scaled) <- diag(covar) / (medias_por_trait ^ 2)   # substitui a diagonal (variâncias / média ^ 2)
  
  # normalizações
  dimor_norm <- dimor / sqrt(sum(dimor ^ 2))
  
  # evolvability
  resposta <- as.matrix(covar) %*% as.vector(dimor_norm)
  evolv <- sum((resposta * dimor_norm) / (size ^ 2))
  #evolv = (diag(t(dimor_norm) %*% covar %*% dimor_norm)) / (dimor_norm ^ 2)
  #lambda <- 1e-6
  #covar_reg <- covar + diag(lambda, nrow(covar))
  #cond = 1 / diag(t(dimor_norm) %*% solve(covar_reg) %*% dimor_norm)
  
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
  integ <- CalcEigenVar(cov2cor(covar))
  
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
    Maior_Lambda = maior_lambda
  )
}

# junta tudo num único data.frame
evolvas <- do.call(rbind, results)
rownames(evolvas) <- NULL

# save results
saveRDS(evolvas, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
#evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")

evolvas <- evolvas[match(tree$tip.label, evolvas$Species), ]
evolvas$genus <- factor(evolvas$Species, levels = tree$tip.label)

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
