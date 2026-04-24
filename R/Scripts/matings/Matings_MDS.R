# MDS Matings
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings")

# load packages
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)){install.packages("tidyr"); library(tidyr)}
if(!require(StatMatch)){install.packages("StatMatch"); library(StatMatch)}
if(!require(FD)){install.packages("FD"); library(FD)}
if(!require(ggrepel)){install.packages("ggrepel"); library(ggrepel)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}
if(!require(reshape2)){install.packages("reshape2"); library(reshape2)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}

# read mating data
dados = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/genus_Matings.csv", sep = ";")
dados = dados[!is.na(dados$DOMINANCE) & !is.na(dados$AGGRESSION), ]

# Classes
#UNIMUNIF(0)  MMMF/UNIMUNIF(1)     MMMF OR UNIMALE/MMMF/UNIMUNIF(2)     UNIMALE/MMMF(3)     UNIMALE(4)
classes = c(unique(dados$SOCIAL_ORGANIZATION))
numeros = c(3, 0, 2, 1, 4, 2)
class_social = data.frame(classes, numeros)

# MONOGAMY(0)  MONOGAMY/POLYGYNADROUS(1)       POLYGYNANDROUS or MONOGAMY/POLYANDROUS or MONOGAMY/POLYGYNY(2)       POLYGYNY/POLYGYNANDROUS (3)         POLYGYNY(4)
classes = c(unique(dados$MATING_SYSTEM))
numeros = c(3, 0, 2, 4, 1, 2, 2)
class_mating = data.frame(classes, numeros)

classes = c(unique(dados$DOMINANCE))
numeros = c(1, 0, 0, 1, 1)
class_dominance = data.frame(classes, numeros)

classes = c(unique(dados$AGGRESSION))
numeros = c(1, 0, 1, 0, 1, 1)
class_aggression = data.frame(classes, numeros)

cols = c("ALL_MALE_GROUPS", "FURTIVE_COPULATION", "INFANTICIDE", "MULTILEVEL_SOCIETY")
dados[cols] <- lapply(dados[cols], function(x) {
  x[is.na(x)] <- "NA"
  as.character(x)
})

classes = c(unique(dados$ALL_MALE_GROUPS))
numeros = c(0, 1)
class_male_groups = data.frame(classes, numeros)

classes = c(unique(dados$FURTIVE_COPULATION))
numeros = c(1, 0)
class_furtive = data.frame(classes, numeros)

classes = c(unique(dados$INFANTICIDE))
numeros = c(1, 0)
class_infanticide = data.frame(classes, numeros)

classes = c(unique(dados$MULTILEVEL_SOCIETY))
numeros = c(0, 1)
class_multilevel = data.frame(classes, numeros)

# transform classes in numbers
list_classes = list(class_social, class_mating, class_dominance, class_aggression,
                    class_male_groups, class_furtive, class_infanticide, class_multilevel)

# Apply classes changes
for(i in 1:2){
  lista = list_classes[[i]]
  coluna = dados[i+3]
  
  for(j in 1:nrow(lista)){
    index = which(coluna == lista[j,1])
    coluna[index,] = lista[j,2]
  }
  dados[i+3] = as.numeric(unlist(coluna))
}

for(i in 3:8){
  lista = list_classes[[i]]
  coluna = dados[i+4]
  
  for(j in 1:nrow(lista)){
    index = which(coluna == lista[j,1])
    coluna[index,] = lista[j,2]
  }
  dados[i+4] = as.numeric(unlist(coluna))
}

#################################
# GOWER DIST + NMDS
#datas = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv", header = TRUE, sep = ";")
genus = dados$GENUS
rownames(dados) = dados$GENUS

dados <- dados %>%
  mutate(full_species = paste(GENUS, SPECIES, sep = "_")) %>%
  filter(full_species %in% species)

# Haplorrhini
fit <- metaMDS(dados[,4:13], distance = "gower", k = 2,
               trymax = 1000, maxit = 1000)
png(
  filename = "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_StressPlot.png",
  width = 4000,        # largura em pixels
  height = 3000,       # altura em pixels
  res = 300            # resolução (dpi)
)
stressplot(fit)
dev.off()

x = as.data.frame(vegan::scores(fit)$sites[,1])
y = as.data.frame(vegan::scores(fit)$sites[,2])

#especies = paste(dados$GENUS, dados$SPECIES, sep = "_")
data = data.frame(x, y, genus, dados$PARVORDER, dados$SOCIAL_ORGANIZATION, dados$MATING_SYSTEM, 
                  dados$PROP_MALES_FEMALES, dados$DOMINANCE, dados$AGGRESSION,dados$ALL_MALE_GROUPS,
                  dados$FURTIVE_COPULATION, dados$INFANTICIDE, dados$MULTILEVEL_SOCIETY, dados$N_PAPERS)

saveRDS(data, "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")

# 1. Extract scores of NMDS
scores_fit = as.data.frame(vegan::scores(fit)$sites) %>%
  mutate(PARVORDER = data$dados.PARVORDER)

# 2. Adjust vectors
fit_env <- envfit(fit, data[, 5:14], permutations = 999)

# 3. Extract coordinate of vectors and apply scale
scale_factor <- 0.5
vectors <- as.data.frame(fit_env$vectors$arrows * sqrt(fit_env$vectors$r) * scale_factor) %>%
  mutate(variable = gsub("^dados\\.", "", rownames(fit_env$vectors$arrows)))  # remove "dados."
novos_nomes <- c("Social Organization", "Mating System", "Proportion M and F", "Dominance", 
                 "Aggression", "All-Male Groups", "Furtive Copulations", "Infanticide", "Multilevel Society", "Number_Papers")  # 9 variáveis
vectors <- vectors %>% mutate(variable = novos_nomes)

# 4. Centroides
centroids <- scores_fit %>% 
  group_by(data$dados.PARVORDER) %>% 
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

parvorders <- unique(data$dados.PARVORDER)
cores <- c("#7570b3", "#1b9e77")
names(cores) <- parvorders

# Plot
mds <- ggplot(scores_fit, aes(x = NMDS1, y = NMDS2, color = PARVORDER)) +
  
  # Elipses com mesmas cores dos pontos
  stat_ellipse(data = scores_fit, aes(x = NMDS1, y = NMDS2, group = PARVORDER, fill = PARVORDER),
              geom = "polygon", alpha = 0.2, color = NA) +

  geom_point(
    data = centroids,
    aes(x = NMDS1, y = NMDS2),  # só as coordenadas
    color = "black",             # cor fixa fora do aes
    shape = 4, size = 5, stroke = 2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_d() +
     
  # Pontos
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2.5, shape = 19, alpha = 0.9, stroke = 2) +
     
  # Vetores
  geom_segment(data = vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  arrow = arrow(length = unit(0.2, "cm")),
  inherit.aes = FALSE, color = "black") +
     
  # Nomes das variáveis
  geom_text_repel(data = vectors,
                  aes(x = NMDS1, y = NMDS2, label = variable),
                      inherit.aes = FALSE, size = 4, max.overlaps = Inf, force = 100,
                      point.padding = 1) +
  
  # Mesmas cores viridis para pontos e elipses
  scale_color_manual(values = cores) +
  scale_fill_manual(values = cores) +
     
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        axis.title = element_text(face = "bold")) +
     
  labs(x = "NMDS1",
       y = "NMDS2",
       color = "Parvorder",
       fill = "Parvorder")

mds       

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.png", plot = mds,
       width = 14,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


tocor = data.frame(scores_fit, dados$SOCIAL_ORGANIZATION, dados$MATING_SYSTEM, dados$PROP_MALES_FEMALES, 
                   dados$DOMINANCE, dados$AGGRESSION, dados$ALL_MALE_GROUPS, dados$FURTIVE_COPULATION,
                   dados$INFANTICIDE, dados$MULTILEVEL_SOCIETY, dados$N_PAPERS)

# selecionar variáveis comportamentais
dados_vars <- tocor[, grep("^dados", names(tocor))]

# matriz com NMDS + dados
nmds_dados <- cbind(tocor[, c("NMDS1", "NMDS2")], dados_vars)

# correlação de Spearman
cor_nmds <- cor(nmds_dados, method = "spearman", use = "pairwise.complete.obs")

# manter apenas NMDS vs dados
cor_nmds[1:2, -c(1:2)]

cor_dados <- cor(dados_vars,
                 method = "spearman",
                 use = "pairwise.complete.obs")

cor_dados

library(Hmisc)

cor_test <- rcorr(as.matrix(nmds_dados), type = "spearman")

cor_test$r   # correlações
cor_test$P   # p-values

