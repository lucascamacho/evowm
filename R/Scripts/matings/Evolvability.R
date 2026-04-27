# Evolvability of sexual dimorphism

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

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

#
matings = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
align = read.table("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv", header = TRUE, sep = ",")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")

# read and plot phylo tree
commom <- intersect(names(vcv), matings$genus)
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
genus <- sub("_.*", "", tree$tip.label)
keep <- genus %in% commom
tree2 <- drop.tip(tree, tree$tip.label[!keep])
genus2 <- sub("_.*", "", tree2$tip.label)
keep_one <- !duplicated(genus2)
tree_genus <- drop.tip(tree2, tree2$tip.label[!keep_one])
tree_genus$tip.label <- sub("_.*", "", tree_genus$tip.label)
#plot(tree_genus)

# cut some vcv
vcv <- vcv[names(vcv) %in% matings$genus]
tree_genus = drop.tip(tree_genus, setdiff(matings$genus, names(vcv)))

all_cov_matrices = PhyloW(tree_genus, vcv)
ancestral = getMRCA(tree_genus, tree_genus$tip.label)

results <- list()  # lista vazia
for(i in seq_along(matings$genus)){
  sp <- matings$genus[i]
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
  # tam_cra <- geomean(medidas$Machos)
  #gm_machos  <- geomean(medidas$Machos)
  #gm_femeas  <- geomean(medidas$Fêmeas)
  #size <- geomean(c(gm_machos, gm_femeas))
  size <- (geomean(medidas$Machos) + geomean(medidas$Fêmeas)) / 2
  #size <- (geomean(log(medidas$Machos)) + geomean(log(medidas$Fêmeas))) / 2
  
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
  
  # evolvability OK
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
  #integ <- CalcR2(cov2cor(covar))
  
  # salva numa linha
  results[[sp]] <- data.frame(
    genus = sp,
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
saveRDS(evolvas, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability.RDS")
#evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")

evolvas <- evolvas[match(tree_genus$tip.label, evolvas$genus), ]
evolvas$genus <- factor(evolvas$genus, levels = tree_genus$tip.label)

p1 = ggplot(evolvas, aes(x = genus, y = Evolvability)) +
     geom_col() +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_discrete(drop = TRUE) +
     xlab(" Genus") +
     ylab("Evolvability SD") +
     coord_flip() +
     theme_classic()

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability_Species.png", plot = p1,
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


# 3 parts Plot 1. No correlation, 2.positive correlation, 3. no correlation
evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Evolvability, y = Size)) +
  #geom_text(aes(x = Evolvability, y = Dimorfism, label = species)) +
  #geom_smooth(aes(x = log(Evolvability), y = Dimorfism), method = "lm") +
  xlab("Evolvability") +
  ylab("Size") +
  #scale_y_continuous(expand = c(0, 0.005)) +
  #scale_x_continuous(expand = c(0, 0.005)) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

summary(lm(evolvas$Dimorfism ~ log(evolvas$Evolvability)))

# Bigger the Avg_Evolvability, bigger the Evolvability
evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Maior_Lambda, y = Evolvability)) +
  xlab("Maior_Lambda") +
  ylab("Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

# Non linear relationship, higher Evolvability, higher SD
evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Evolvability, y = Dimorfism)) +
  xlab("Evolvability") +
  ylab("SD") +
  #xlim(0, 0.2) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)# tamanho do título eixo y
  )

evolv_plot

# Linear, Higher Integration, Higher SD
integ_sd_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Integration, y = Dimorfism)) +
  geom_smooth(aes(x = Integration, y = Dimorfism), method = "lm") +
  xlab("Integration") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

integ_sd_plot

summary(lm(evolvas$Dimorfism ~ evolvas$Integration))

# Linear, High Integration, High Evolvability
integ_evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Integration, y = Evolvability)) +
  xlab("Integration") +
  ylab("Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

integ_evolv_plot

# With species names: High Evolvability, High SD
evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Evolvability, y = Dimorfism, label = species)) +
  xlab("Evolvability") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

# Non-linear: High Conditional Evolvability, High SD
evolv_plot <- ggplot(data = evolvas) +
  #geom_point(aes(x = Conditional_Evolvability, y = Dimorfism)) +
  geom_text(aes(x = Conditional_Evolvability, y = Dimorfism, label = species)) +
  xlab("Conditional_Evolvability") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

# Higher Conditional Evolvability, Higher Evolvability But there is different cases
evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Conditional_Evolvability, y = Evolvability, label = species)) +
  xlab("Conditional_Evolvability") +
  ylab("Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

# Non-Linear: Higher Integration, Lower Conditional Evolvability
evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Integration_A, y = Conditional_Evolvability, label = species)) +
  xlab("Integration_A") +
  ylab("Conditional_Evolvability") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

# No correlation between Size and Integration
evolv_plot <- ggplot(data = evolvas) +
  geom_text(aes(x = Size, y = Integration_A, label = species)) +
  xlab("Size") +
  ylab("Integration_A") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

evolv_plot

# Linear, Higher Integration, Higher SD
integ_sd_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Integration_A, y = Dimorfism)) +
#  geom_smooth(aes(x = Integration_A, y = Dimorfism), method = "lm") +
  xlab("Integration_A") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

integ_sd_plot

# Sexual Dimorphism and Size
integ_sd_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Size, y = Dimorfism)) +
  #  geom_smooth(aes(x = Integration_A, y = Dimorfism), method = "lm") +
  xlab("Size") +
  ylab("SD") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12), # tamanho dos números no eixo x
    axis.text.y  = element_text(size = 12),                                     # tamanho dos números no eixo y
    axis.title.x = element_text(size = 14),                                     # tamanho do título eixo x
    axis.title.y = element_text(size = 14)                                      # tamanho do título eixo y
  )

integ_sd_plot

#
# Non-linear Higher alignments between Pmax and Dimorphism, Higher the magnitude of sexual dimorphism
dimor_data = data.frame(evolvas$species, evolvas$Standart_Evolvability, evolvas$Dimorfism, 
                        evolvas$Conditional_Evolvability, evolvas$Integration_A, evolvas$Evolvability,
                        align$align_1)

dimor_plot = ggplot(data = dimor_data) +
  geom_point(aes(x = align.align_1, y = evolvas.Dimorfism)) +
  geom_smooth(aes(x = align.align_1, y = evolvas.Dimorfism), method = "lm", se = FALSE, formula = y ~ x) +
  xlab("Pmax x Sexual Dimorphism Vector") +
  ylab("Standardized Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),  # Tamanho das legendas dos eixos
    axis.text = element_text(size = 12)    # Tamanho dos números dos eixos
  )

dimor_plot

#ggsave(dimor_plot_2, filename = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Alinhamento_Dimorphism_Evolvability.pdf", dpi = 600,
#       width = 35, height = 20, units = "cm",  bg = "transparent")
