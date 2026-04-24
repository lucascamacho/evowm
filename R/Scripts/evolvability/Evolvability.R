# Evolvability of sexual dimorphism
# Several measures of evolvability
# in log, removing PTTSP
# and removing the first eingenvalue of the variance matrix
# Also, using the pooled matrix of all species (ancestor species)

setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability")

library(ggplot2)
library(openxlsx)
library(evolqg)
library(evolvability)
library(ape)
library(ggpmisc)

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# read V/CV matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
names(vcv)  = gsub(".csv", replacement= "", temp)

#
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")

#
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
species = names(vcv)
species = append(species, "Homo_sapiens")
tree = drop.tip(tree, setdiff(tree$tip.label, species))

#
all_cov_matrices = PhyloW(tree, vcv)
ancestral = getMRCA(tree, tree$tip.label)
covar = all_cov_matrices[[ancestral]]

#
eig = eigen(covar)
D = diag(eig$values)
V = eig$vectors
D2 = D
D2[1,1] = 0
covar = V %*% D2 %*% t(V)

# colocar humanos
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")

# check names and remove doubts
msrs$SEX[which(msrs$SEX == "?female")] = "female"
msrs$SEX[which(msrs$SEX == "?male")] = "male"

# remove uncertain sex
msrs$SEX[msrs$SEX == "0"] = NA
msrs$SEX[msrs$SEX == ""] = NA
msrs$SEX[msrs$SEX == "sexo"] = NA

msrs = msrs[complete.cases(msrs$SEX), ]

species_subset = msrs[which(msrs$GENUS == "Homo"), ]
species_subset = species_subset[which(species_subset$SPECIES == "sapiens"), ]

# separate M and F
sub_sexes_m = species_subset[which(species_subset$SEX == "male"), ]
sub_sexes_f = species_subset[which(species_subset$SEX == "female"), ]

# geometric means of M and F
Machos = apply(log(sub_sexes_m[,c(49:67, 69:87)]), 2, geomean)
Fêmeas = apply(log(sub_sexes_f[,c(49:67, 69:87)]), 2, geomean)

medias$ByTrait_Averages$Homo_sapiens <- list(
  Machos = Machos,
  Fêmeas = Fêmeas
)

#
results <- list()  # lista vazia
for(i in seq_along(species)){
  sp <- species[i]
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
  
  # 
  size <- (geomean(medidas$Machos) + geomean(medidas$Fêmeas)) / 2

  #
  dimor <- medidas$Machos - medidas$Fêmeas

  #
  covar <- vcv[[sp]]
  covar <- as.matrix(covar)

  # remove 1rst eingenvector
  eig = eigen(covar)
  D = diag(eig$values)
  V = eig$vectors
  D2 = D
  D2[1,1] = 0
  covar = V %*% D2 %*% t(V)
  
  # normalizações
  dimor_norm <- dimor / sqrt(sum(dimor ^ 2))
  
  # evolvability OK
  evolv <- as.numeric(t(dimor_norm) %*% covar %*% dimor_norm)

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
  integ <- CalcEigenVar(as.matrix(covar))
  
  # salva numa linha
  results[[sp]] <- data.frame(
    species = sp,
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

evolvas <- evolvas[match(tree$tip.label, evolvas$species), ]
evolvas$species <- factor(evolvas$species, levels = tree$tip.label)

# save results
saveRDS(evolvas, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability.RDS")
#evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")

p1 = ggplot(evolvas, aes(x = species, y = Evolvability)) +
     geom_col() +
     scale_y_continuous(expand = c(0, 0)) +
     scale_x_discrete(drop = TRUE) +
     xlab(" Species") +
     ylab("Evolvability SD") +
     coord_flip() +
     theme_classic()

p1
#ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Species.png", plot = p1,
#       width = 12,    # largura em inches
#       height = 8,   # altura em inches
#       dpi = 200)    # resolução

evolv_plot = ggplot(data = evolvas) +
  geom_point(aes(x = Dimorfism, y = Evolvability), size = 2) +
  stat_poly_eq(
    aes(x = Dimorfism, y = Evolvability, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 5,
    label.x = "right"
  ) +
  geom_smooth(aes(x = Dimorfism, y = Evolvability), method = "lm", se = TRUE) +
  xlab("Sexual Dimorfism") +
  ylab("Evolvability Estimates of Sexual Dimorphism") +
  theme_classic() +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )

evolv_plot

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Dimorphism.png", plot = evolv_plot,
       width = 12,
       height = 7,
       dpi = 300
)

