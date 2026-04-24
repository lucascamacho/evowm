# vector correlation between PCs and isometric modules vectors

# load packages and functions
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
library(dplyr)
library(tidyr)

geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# load data to correlation
medidas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
matings = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")

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
comum <- intersect(names(vcv), medidas$Genus)

# create modules and isometric vectors 
mod_names = c("ISPM", "ISNSL", "ISPNS", "PMZS", "PMZI", "PMMT", "NSLNA", "NSLZS", 
              "NSLZI", "NABR", "NAFM", "NAPNS", "BRPT", "BRAPET", "PTFM", "PTAPET",
              "PTBA", "PTEAM", "PTZYGO", "PTTSP", "FMZS", "FMMT", "ZSZI", "ZIMT",
              "ZIZYGO", "ZITSP", "MTPNS", "PNSAPET", "APETBA", "APETTS", "BAEAM", 
              "EAMZYGO", "ZYGOTSP", "LDAS", "BRLD", "OPILD", "PTAS", "JPAS", "BAOPI")

isometric = rep(0.160128154, 39)
modules = c(1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, 
            -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, -1, -1)
names(modules) = mod_names

cor_iso_1 = vector()
cor_iso_2 = vector()
cor_iso_3 = vector()
cor_iso_4 = vector()
cor_iso_5 = vector()
cor_iso_6 = vector()
cor_iso_7 = vector()
cor_iso_8 = vector()

cor_modules_1 = vector()
cor_modules_2 = vector()
cor_modules_3 = vector()
cor_modules_4 = vector()
cor_modules_5 = vector()
cor_modules_6 = vector()
cor_modules_7 = vector()
cor_modules_8 = vector()

for(i in 1:length(comum)){
  # chose species
  sp = comum[i]
  
  # pc
  pc_1 = medidas$PCs[[sp]][1]
  pc_2 = medidas$PCs[[sp]][2]
  pc_3 = medidas$PCs[[sp]][3]
  pc_4 = medidas$PCs[[sp]][4]
  pc_5 = medidas$PCs[[sp]][5]
  pc_6 = medidas$PCs[[sp]][6]
  pc_7 = medidas$PCs[[sp]][7]
  pc_8 = medidas$PCs[[sp]][8]
  
  # cors iso
  cor_iso_1[i] = abs(corVector(unlist(pc_1), isometric))
  cor_iso_2[i] = abs(corVector(unlist(pc_2), isometric))
  cor_iso_3[i] = abs(corVector(unlist(pc_3), isometric))
  cor_iso_4[i] = abs(corVector(unlist(pc_4), isometric))
  cor_iso_5[i] = abs(corVector(unlist(pc_5), isometric))
  cor_iso_6[i] = abs(corVector(unlist(pc_6), isometric))
  cor_iso_7[i] = abs(corVector(unlist(pc_7), isometric))
  cor_iso_8[i] = abs(corVector(unlist(pc_8), isometric))
  
  # cors modules
  cor_modules_1[i] = abs(corVector(unlist(pc_1), modules))
  cor_modules_2[i] = abs(corVector(unlist(pc_2), modules))
  cor_modules_3[i] = abs(corVector(unlist(pc_3), modules))
  cor_modules_4[i] = abs(corVector(unlist(pc_4), modules))
  cor_modules_5[i] = abs(corVector(unlist(pc_5), modules))
  cor_modules_6[i] = abs(corVector(unlist(pc_6), modules))
  cor_modules_7[i] = abs(corVector(unlist(pc_7), modules))
  cor_modules_8[i] = abs(corVector(unlist(pc_8), modules))
}

atuais = data.frame(comum, cor_iso_1, cor_iso_2, cor_iso_3, cor_iso_4, cor_iso_5, cor_iso_6, cor_iso_7, cor_iso_8,
                    cor_modules_1, cor_modules_2, cor_modules_3, cor_modules_4, cor_modules_5, cor_modules_6, cor_modules_7, cor_modules_8)

# ancestrals
medidas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")

cor_iso_1 = vector()
cor_iso_2 = vector()
cor_iso_3 = vector()
cor_iso_4 = vector()
cor_iso_5 = vector()
cor_iso_6 = vector()
cor_iso_7 = vector()
cor_iso_8 = vector()

cor_modules_1 = vector()
cor_modules_2 = vector()
cor_modules_3 = vector()
cor_modules_4 = vector()
cor_modules_5 = vector()
cor_modules_6 = vector()
cor_modules_7 = vector()
cor_modules_8 = vector()

for(i in 1:length(medidas$Ancestrals)){
  # chose species
  sp = medidas$Ancestrals[[i]]
  index = which(names(medidas$VCV) == sp)
  
  # pc
  pc_1 = medidas$PCs[[as.character(sp)]][1]
  pc_2 = medidas$PCs[[as.character(sp)]][2]
  pc_3 = medidas$PCs[[as.character(sp)]][3]
  pc_4 = medidas$PCs[[as.character(sp)]][4]
  pc_5 = medidas$PCs[[as.character(sp)]][5]
  pc_6 = medidas$PCs[[as.character(sp)]][6]
  pc_7 = medidas$PCs[[as.character(sp)]][7]
  pc_8 = medidas$PCs[[as.character(sp)]][8]
  
  # cors iso
  cor_iso_1[i] = abs(corVector(unlist(pc_1), isometric))
  cor_iso_2[i] = abs(corVector(unlist(pc_2), isometric))
  cor_iso_3[i] = abs(corVector(unlist(pc_3), isometric))
  cor_iso_4[i] = abs(corVector(unlist(pc_4), isometric))
  cor_iso_5[i] = abs(corVector(unlist(pc_5), isometric))
  cor_iso_6[i] = abs(corVector(unlist(pc_6), isometric))
  cor_iso_7[i] = abs(corVector(unlist(pc_7), isometric))
  cor_iso_8[i] = abs(corVector(unlist(pc_8), isometric))
  
  # cors modules
  cor_modules_1[i] = abs(corVector(unlist(pc_1), modules))
  cor_modules_2[i] = abs(corVector(unlist(pc_2), modules))
  cor_modules_3[i] = abs(corVector(unlist(pc_3), modules))
  cor_modules_4[i] = abs(corVector(unlist(pc_4), modules))
  cor_modules_5[i] = abs(corVector(unlist(pc_5), modules))
  cor_modules_6[i] = abs(corVector(unlist(pc_6), modules))
  cor_modules_7[i] = abs(corVector(unlist(pc_7), modules))
  cor_modules_8[i] = abs(corVector(unlist(pc_8), modules))
}

ancestrais = data.frame(medidas$Ancestrals, cor_iso_1, cor_iso_2, cor_iso_3, cor_iso_4, cor_iso_5, cor_iso_6, cor_iso_7, cor_iso_8,
                              cor_modules_1, cor_modules_2, cor_modules_3, cor_modules_4, cor_modules_5, cor_modules_6, cor_modules_7, cor_modules_8)

# Padronizar o nome da coluna
colnames(atuais)[1] <- "medidas.Ancestrals"

# Unir os data frames
novo_df <- rbind(atuais, ancestrais)

names(novo_df) = c("Species", "cor_iso_1", "cor_iso_2", "cor_iso_3", "cor_iso_4", "cor_iso_5", "cor_iso_6", "cor_iso_7", "cor_iso_8",
                            "cor_modules_1", "cor_modules_2", "cor_modules_3", "cor_modules_4", "cor_modules_5", "cor_modules_6", "cor_modules_7", "cor_modules_8")

write.csv(novo_df, "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_Iso_Mod.csv", row.names = FALSE)

#
data = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_iso_mod.csv")

cols <- 2:5
titles <- paste("Isometric and PC", 1:4)

# Coloca os dados em formato longo
df_long <- data %>%
  select(all_of(cols)) %>%
  setNames(titles) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Component",
    values_to = "Correlation"
  )

# Cria o gráfico
# Cria o gráfico com histograma + densidade
p <- ggplot(df_long, aes(x = Correlation)) +
  geom_histogram(
    aes(y = ..density..),
    bins = 20,
    fill = "#74c476",    # verde claro (agradável)
    color = "white",
    alpha = 0.8
  ) +
  geom_density(
    color = "#238b45",   # verde escuro
    linewidth = 1.2,
    alpha = 0.8
  ) +
  facet_wrap(~ Component, scales = "free", ncol = 4) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Correlation",
    y = "Density",
    title = "Alignment Between Principal Components and Isometric Vector"
  )
p

# Salva o gráfico em alta resolução
ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Hist_All_cor_PCS_Isometric.png",
  plot = p,
  width = 12,
  height = 7,
  dpi = 300
)

