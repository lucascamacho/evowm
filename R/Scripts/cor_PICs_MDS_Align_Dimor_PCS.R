# set WD and functions
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

if(!require(ape)){install.packages("ape"); library(ape)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)){install.packages("tidyr"); library(tidyr)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}


# load data
matings = read.table("~/Dropbox/Doc/Data/wos_mating_systems/Matings.csv", sep = ",", header = TRUE)
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS") 
dimor = read.table("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv", header = TRUE, sep = ",", row.names = 1)

# Catarrhini
species = dimor$matings.especies[1:38]

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
tree = drop.tip(tree, setdiff(tree$tip.label, species))

# cut mds and dimor to only the catarrhini
mds_filtrado = mds[mds$especies %in% species, ]
dimor_filtrado = dimor[dimor$matings.especies %in% species, ]

# Vetores de colunas a testar
cols_dimor <- 4:11
cols_mds <- 1:2

# Cria data frame vazio
co <- data.frame()

# Loop sobre combinações de colunas
for (mds_col in cols_mds) {
  for (dimor_col in cols_dimor) {
    r <- cor.test(dimor_filtrado[[dimor_col]], mds_filtrado[[mds_col]])
    co <- rbind(co, data.frame(
      dimor_col = names(dimor_filtrado)[dimor_col],
      mds_col = mds_col,
      estimate = r$estimate,
      p.value = r$p.value
    ))
  }
}

# calculate fisher correlations
co$z_fisher = 0.5 * log((1 + co$estimate) / (1 - co$estimate))

# p value for fisher
n = nrow(mds_filtrado)
z_stat = co$z_fisher * sqrt(n - 3)
co$p_z_fisher = 2 * (1 - pnorm(abs(z_stat)))


# PIC para NMDS1 (eixos x)
pics_nmds1 = vector()
p_values_pics_nmds1 = vector()
for(i in 1:8){
  mds_1 = setNames(mds_filtrado$scores.fit....1., mds_filtrado$especies)
  dimor_pc = setNames(dimor_filtrado[,i+3], dimor_filtrado$Species)
  
  pic_mds = pic(mds_1, tree)
  pic_dimor = pic(dimor_pc, tree)
  
  r = cor(pic_mds, pic_dimor)
  z = 0.5 * log((1 + r) / (1 - r))
  
  pics_nmds1[i] = z
  
  zse = 1 / sqrt(length(pic_mds) - 3)
  z_stat = z / zse
  p_values_pics_nmds1[i] = 2 * (1 - pnorm(abs(z_stat)))}

# PIC para NMDS2 (eixos y)
pics_nmds2 = vector()
p_values_pics_nmds2 = vector()
for(i in 1:8){
  mds_2 = setNames(mds_filtrado$scores.fit....2., mds_filtrado$especies)
  dimor_pc = setNames(dimor_filtrado[,i+3], dimor_filtrado$Species)
  
  pic_mds = pic(mds_2, tree)
  pic_dimor = pic(dimor_pc, tree)
  
  r = cor(pic_mds, pic_dimor)
  z = 0.5 * log((1 + r) / (1 - r))
  
  pics_nmds2[i] = z
  
  zse = 1 / sqrt(length(pic_mds) - 3)
  z_stat = z / zse
  p_values_pics_nmds2[i] = 2 * (1 - pnorm(abs(z_stat)))
}

co = cbind(co, pics_nmds1, p_values_pics_nmds1, pics_nmds2, p_values_pics_nmds2)

# getting all together
colnames(co)[1:4] = c("Align", "NMDS", "simple_r", "r_p_value")
row.names(co) = NULL
co_c = co

# Platyrrhini
species = dimor$matings.especies[39:62]

# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
tree = drop.tip(tree, setdiff(tree$tip.label, species))

# cut mds and dimor to only the catarrhini
mds_filtrado = mds[mds$especies %in% species, ]
dimor_filtrado = dimor[dimor$matings.especies %in% species, ]

# Vetores de colunas a testar
cols_dimor <- 4:11
cols_mds <- 1:2

# Cria data frame vazio
co <- data.frame()

# Loop sobre combinações de colunas
for (mds_col in cols_mds) {
  for (dimor_col in cols_dimor) {
    r <- cor.test(dimor_filtrado[[dimor_col]], mds_filtrado[[mds_col]])
    co <- rbind(co, data.frame(
      dimor_col = names(dimor_filtrado)[dimor_col],
      mds_col = mds_col,
      estimate = r$estimate,
      p.value = r$p.value
    ))
  }
}

# calculate fisher correlations
co$z_fisher = 0.5 * log((1 + co$estimate) / (1 - co$estimate))

# p value for fisher
n = nrow(mds_filtrado)
z_stat = co$z_fisher * sqrt(n - 3)
co$p_z_fisher = 2 * (1 - pnorm(abs(z_stat)))


# PIC para NMDS1 (eixos x)
pics_nmds1 = vector()
p_values_pics_nmds1 = vector()
for(i in 1:8){
  mds_1 = setNames(mds_filtrado$scores.fit....1., mds_filtrado$especies)
  dimor_pc = setNames(dimor_filtrado[,i+3], dimor_filtrado$Species)
  
  pic_mds = pic(mds_1, tree)
  pic_dimor = pic(dimor_pc, tree)
  
  r = cor(pic_mds, pic_dimor)
  z = 0.5 * log((1 + r) / (1 - r))
  
  pics_nmds1[i] = z
  
  zse = 1 / sqrt(length(pic_mds) - 3)
  z_stat = z / zse
  p_values_pics_nmds1[i] = 2 * (1 - pnorm(abs(z_stat)))
}

# PIC para NMDS2 (eixos y)
pics_nmds2 = vector()
p_values_pics_nmds2 = vector()
for(i in 1:8){
  mds_2 = setNames(mds_filtrado$scores.fit....2., mds_filtrado$especies)
  dimor_pc = setNames(dimor_filtrado[,i+3], dimor_filtrado$Species)
  
  pic_mds = pic(mds_2, tree)
  pic_dimor = pic(dimor_pc, tree)
  
  r = cor(pic_mds, pic_dimor)
  z = 0.5 * log((1 + r) / (1 - r))
  
  pics_nmds2[i] = z
  
  zse = 1 / sqrt(length(pic_mds) - 3)
  z_stat = z / zse
  p_values_pics_nmds2[i] = 2 * (1 - pnorm(abs(z_stat)))
}

co = cbind(co, pics_nmds1, p_values_pics_nmds1, pics_nmds2, p_values_pics_nmds2)

# getting all together
colnames(co)[1:4] = c("Align", "NMDS", "simple_r", "r_p_value")
row.names(co) = NULL
co_p = co

# save
saveRDS(list(co_c, co_p), file = "~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PICs_MDS_Align_Dimor_PCS.RDS")

# Plots
# Catarrhini
co_c$pics_nmds1_r <- tanh(co_c$pics_nmds1)
co_c$pics_nmds2_r <- tanh(co_c$pics_nmds2)

paraplot_c = data.frame(
  Align = co_c$Align,
  NMDS_orig = co_c$NMDS,  # renomeia a coluna que existia
  pics_nmds1 = co_c$pics_nmds1_r,   
  p_values_nmds1 = co_c$p_values_pics_nmds1,
  pics_nmds2 = co_c$pics_nmds2_r,   
  p_values_nmds2 = co_c$p_values_pics_nmds2
)

# Dados já transformados para formato longo
co_long <- paraplot_c %>%
  pivot_longer(
    cols = c(pics_nmds1, pics_nmds2),
    names_to = "NMDS",      # agora não entra em conflito
    values_to = "estimate"
  ) %>%
  mutate(
    p_value = ifelse(NMDS == "pics_nmds1", p_values_nmds1, p_values_nmds2),
    NMDS = ifelse(NMDS == "pics_nmds1", "NMDS1", "NMDS2"),
    label = sprintf("%.2f", estimate),
    fontface = case_when(
      p_value < 0.01 ~ "bold",
      p_value < 0.05 ~ "italic",
      TRUE ~ "plain"
    ),
    Align = paste0("Align_SD_PC", as.numeric(gsub("align_", "", Align)))
  )


# Plot com ajustes estéticos
c = ggplot(co_long, aes(x = Align, y = NMDS, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label, fontface = fontface), size = 5) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Correlation"
  ) +
  labs(x = NULL, y = NULL, title = "Catarrhini") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "none"
  )


# Platyrrhini
co_p$pics_nmds1_r <- tanh(co_p$pics_nmds1)
co_p$pics_nmds2_r <- tanh(co_p$pics_nmds2)

paraplot_p = data.frame(
  Align = co_p$Align,
  NMDS_orig = co_p$NMDS,  # renomeia a coluna que existia
  pics_nmds1 = co_p$pics_nmds1_r,   
  p_values_nmds1 = co_p$p_values_pics_nmds1,
  pics_nmds2 = co_p$pics_nmds2_r,   
  p_values_nmds2 = co_p$p_values_pics_nmds2
)

co_long <- paraplot_p %>%
  pivot_longer(
    cols = c(pics_nmds1, pics_nmds2),
    names_to = "NMDS",      # agora não entra em conflito
    values_to = "estimate"
  ) %>%
  mutate(
    p_value = ifelse(NMDS == "pics_nmds1", p_values_nmds1, p_values_nmds2),
    NMDS = ifelse(NMDS == "pics_nmds1", "NMDS1", "NMDS2"),
    label = sprintf("%.2f", estimate),
    fontface = case_when(
      p_value < 0.01 ~ "bold",
      p_value < 0.05 ~ "italic",
      TRUE ~ "plain"
    ),
    Align = paste0("Align_SD_PC", as.numeric(gsub("align_", "", Align)))
  )


# Plot com ajustes estéticos
p = ggplot(co_long, aes(x = Align, y = NMDS, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label, fontface = fontface), size = 5) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Correlation"
  ) +
  labs(x = NULL, y = NULL, title = "Platyrrhini") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

fig_final = c + p
fig_final

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PICs_MDS_Align_Dimor.png", plot = fig_final,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução
