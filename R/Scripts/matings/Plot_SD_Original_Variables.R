setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

library(ape)
library(ggplot2)
library(dplyr)
library(ggpmisc)
library(geiger)
library(purrr)
library(gt)

# datas
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")

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

#Load and plot the phylogeny
# read and plot phylo tree
commom <- intersect(names(vcv), mds$genus)
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
genus <- sub("_.*", "", tree$tip.label)
keep <- genus %in% commom
tree2 <- drop.tip(tree, tree$tip.label[!keep])
tree <- tree2
# extrair generos
genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
unique_genera <- unique(genera)
representatives <- c()
for(g in unique_genera){
  
  spp <- tree$tip.label[genera == g]
  
  if(length(spp) == 1){
    
    representatives[g] <- spp
    
  } else {
    
    # pegar comprimento do ramo terminal
    edges <- match(spp, tree$tip.label)
    terminal_lengths <- tree$edge.length[match(edges, tree$edge[,2])]
    
    # escolher o menor (mais recente)
    representatives[g] <- spp[which.max(terminal_lengths)]
  }
}

# remover todas as outras espécies
tree_genus <- drop.tip(tree, setdiff(tree$tip.label, representatives))

# renomear para os gêneros
tree_genus$tip.label <- names(representatives)

df = data.frame(mds$genus, mds$dados.PARVORDER, mds$dados.SOCIAL_ORGANIZATION, mds$dados.MATING_SYSTEM,
                mds$dados.PROP_MALES_FEMALES, mds$dados.DOMINANCE, mds$dados.AGGRESSION, mds$dados.ALL_MALE_GROUPS,
                mds$dados.FURTIVE_COPULATION, mds$dados.INFANTICIDE, mds$dados.MULTILEVEL_SOCIETY, mds$dados.N_PAPERS,
                evolvas$Size, align$normas, align$normas, evolvas$Integration, align$dimor)


mycols <- c("#1b9e77", "#7570b3")

df_c = df[which(df$mds.dados.PARVORDER == "Catarrhini"),]
p_c = summary(lm(align.normas ~ evolvas.Integration, data = df_c))$r.squared

df_p = df[which(df$mds.dados.PARVORDER == "Platyrrhini"),]
p_p = summary(lm(align.normas ~ evolvas.Integration, data = df_p))$r.squared

p1 = ggplot(df, aes(x = evolvas.Integration,
               y = align$normas,
               color = mds.dados.PARVORDER, fill = mds.dados.PARVORDER)) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.4,
              linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6, col = "red") +
  theme_bw(base_size = 14) +
  theme_classic(base_size = 14) +
  stat_poly_eq(
    aes(label = paste(..rr.label..), color = mds.dados.PARVORDER),
    formula = y ~ x,
    parse = TRUE,
    size = 6,
    label.x = "left",
    label.y = "bottom"
  ) +
  labs(
    x = "Integration Index",
    y = "Sexual Dimorphism",
    color = "Parvorder"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14)
  ) +
  scale_color_manual(values = mycols) +
  scale_fill_manual(values = mycols) +
  guides(fill = "none")

p1

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_Dimorphism_Extant_integrationIndex.png", plot = p1,
#        width = 12,    # largura em inches
#       height = 8,   # altura em inches
#        dpi = 200)    # resolução


# ANOVA in SOCIAL ORG AND SD
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  social = df$mds.dados.SOCIAL_ORGANIZATION,
  dimorfismo = df$align.normas
)

labels_social = c(1,2,3,4,5)

#labels_social = c("UNIMUNIF", "MMMF/UNIMUNIF","MMMF or UNIMALE/UNIMUNIF or UNIMALE/MMMF/UNIMUNIF or UNIMUNIF/UNIFEMALE",
#                  "MMMF/UNIFEMALE or MMMF/UNIMALE", "UNIMALE or UNIFEMALE")

df_anova$social <- factor(df_anova$social,
                          labels = labels_social)

df_anova$social <- factor(df_anova$social,
                          levels = sort(unique(df_anova$social)),
                          labels = labels_social)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
tree_c = drop.tip(tree_genus, setdiff(tree_genus$tip.label, df_anova_c$genus))
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$social
names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

c_social = aov.phylo(dimor ~ soc, tree_c)
summary(c_social)

df_anova_c = data.frame(dimor, soc)
c = ggplot(df_anova_c, aes(x = soc, y = dimor)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Social Organization",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_SOCIAL_ORG(0.08).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
tree_p = drop.tip(tree_genus, setdiff(tree_genus$tip.label, df_anova_p$genus))
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$social

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

p_social = aov.phylo(dimor ~ soc, tree_p)
summary(p_social)

p = ggplot(df_anova_p, aes(x = social, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Social Organization",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_SOCIAL_OR(0.0038)G.png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

# ANOVA in MATING SYSTEMS AND SD
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  matings = df$mds.dados.MATING_SYSTEM,
  dimorfismo = df$align.normas
)

labels_matings = c(0,1,2,3,4)

#labels_matings = c("MONOGAMY", "POLYGYNANDROUS or  MONOGAMY/POLYGYNY or MONOGAMY/POLYANDROUS", "POLYGYNY/POLYGYNANDROUS",
#                   "POLYGYNY or POLYANDROUS")

df_anova$matings <- factor(df_anova$matings,
                          labels = labels_matings)

df_anova$matings <- factor(df_anova$matings,
                          levels = sort(unique(df_anova$matings)),
                          labels = labels_matings)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$matings

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

c_matings = aov.phylo(dimor ~ soc, tree_c)
summary(c_matings)

c = ggplot(df_anova_c, aes(x = matings, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Mating Systems",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c


# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_MATING_SYSTEM(0.10).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$matings

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

p_matings = aov.phylo(dimor ~ soc, tree_p)
summary(p_matings)

p = ggplot(df_anova_p, aes(x = matings, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Mating Systems",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_MATING_SYSTEM(0.01).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

# LM in PROP M/F AND SD
df_lm <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  prop = df$mds.dados.PROP_MALES_FEMALES,
  dimorfismo = df$align.normas
)

df_lm_c = df_lm[which(df_lm$parvorder == "Catarrhini"),]
c_lm_model <- lm(dimorfismo ~ prop, data = df_lm_c)
r2 <- summary(c_lm_model)$r.squared %>% round(3)
summary(c_lm_model)

c = ggplot(df_lm_c, aes(x = prop, y = dimorfismo)) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2, color = "#1b9e77") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,  # canto superior direito
           label = paste0("R² = ", r2),
           hjust = 1.1, vjust = 1.5,
           size = 5,
           color = "#1b9e77") +
  labs(
    x = "Proportion M/F",
    y = "Sexual Dimorphism"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_blank()
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_PROP_(0.333).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

df_lm_p = df_lm[which(df_lm$parvorder == "Platyrrhini"),]
p_lm_model <- lm(dimorfismo ~ prop, data = df_lm_p)
r2 <- summary(p_lm_model)$r.squared %>% round(3)
summary(p_lm_model)

p = ggplot(df_lm_p, aes(x = prop, y = dimorfismo)) +
  geom_point(size = 4, alpha = 0.9, stroke = 2, color = "#7570b3") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
  annotate("text",
           x = Inf, y = Inf,  # canto superior direito
           label = paste0("R² = ", r2),
           hjust = 1.1, vjust = 1.5,
           size = 5,
           color = "#7570b3") +
  labs(
    x = "Proportion M/F",
    y = "Sexual Dimorphism"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_PROP(0.02).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução


# ANOVA in DOMINANCE AND SD
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  dominance = df$mds.dados.DOMINANCE,
  dimorfismo = df$align.normas
)

labels_dominance = c("NOT_APPLY or NOT_OBSERVED", "MALE, MALE_FEMALE, FEMALE")

df_anova$dominance <- factor(df_anova$dominance,
                           labels = labels_dominance)

df_anova$dominance <- factor(df_anova$dominance,
                           levels = sort(unique(df_anova$dominance)),
                           labels = labels_dominance)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$dominance

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
c_dominance <- aov.phylo(dimor ~ soc, tree_c)
summary(c_dominance)

c = ggplot(df_anova_c, aes(x = dominance, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Dominance",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_DOMINANCE_(0.069).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$dominance

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
p_dominance <- aov.phylo(dimor ~ soc, tree_p)
summary(p_dominance)

p = ggplot(df_anova_p, aes(x = dominance, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Dominance",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_DOMINANCE(0.47).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

# ANOVA in AGGRESSION AND SD
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  aggression = df$mds.dados.AGGRESSION,
  dimorfismo = df$align.normas
)

labels_aggression = c("NOT_OBSERVED", "OBSERVED")

df_anova$aggression <- factor(df_anova$aggression,
                             labels = labels_aggression)

df_anova$aggression <- factor(df_anova$aggression,
                             levels = sort(unique(df_anova$aggression)),
                             labels = labels_aggression)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$aggression

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
c_aggression <- aov.phylo(dimor ~ soc, tree_c)
summary(c_aggression)

c = ggplot(df_anova_c, aes(x = aggression, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Aggression",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_AGGRESSION(0.069).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução


df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$aggression

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
p_aggression <- aov.phylo(dimor ~ soc, tree_p)
summary(p_aggression)

p = ggplot(df_anova_p, aes(x = aggression, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Sexual Dimorphism",
       y = "Aggression") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_AGGRESSION(0.62).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

# ANOVA in ALL MALE GROUPS AND SD
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  all_male = df$mds.dados.ALL_MALE_GROUPS,
  dimorfismo = df$align.normas
)

labels_all_male = c("ABSENCE", "PRESENCE")

df_anova$all_male <- factor(df_anova$all_male,
                              labels = labels_all_male)

df_anova$all_male <- factor(df_anova$all_male,
                              levels = sort(unique(df_anova$all_male)),
                              labels = labels_all_male)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$all_male

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
c_all_male <- aov.phylo(dimor ~ soc, tree_c)
summary(c_all_male)

c = ggplot(df_anova_c, aes(x = all_male, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of All Male Groups",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_ALLMALE(0.032).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução


df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$all_male

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
p_all_male <- aov.phylo(dimor ~ soc, tree_p)
summary(p_all_male)

p = ggplot(df_anova_p, aes(x = all_male, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of All Male Groups",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )
p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_ALLMALE(0.80).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

# ANOVA in FURTIVE COPULATION
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  furtive = df$mds.dados.FURTIVE_COPULATION,
  dimorfismo = df$align.normas
)

labels_furtive = c("ABSENCE", "PRESENCE")

df_anova$furtive <- factor(df_anova$furtive,
                            labels = labels_furtive)

df_anova$furtive <- factor(df_anova$furtive,
                            levels = sort(unique(df_anova$furtive)),
                            labels = labels_furtive)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$furtive

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
c_furtive <- aov.phylo(dimor ~ soc, tree_c)
summary(c_furtive)

c = ggplot(df_anova_c, aes(x = furtive, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Furtive Copulations",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_FURTIVE(0.047).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução


df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$furtive

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
p_furtive <- aov.phylo(dimor ~ soc, tree_p)
summary(p_furtive)

p = ggplot(df_anova_p, aes(x = furtive, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Furtive Copulation",
       y = "Sexual Dimorphism") +
  #         label = round(p_value, 3), size = 5) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_FURTIVE(0.058).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

# ANOVA IN INFANTICIDE AND SD
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  infanticide = df$mds.dados.INFANTICIDE,
  dimorfismo = df$align.normas
)

labels_infanticide = c("ABSENCE", "PRESENCE")

df_anova$infanticide <- factor(df_anova$infanticide,
                           labels = labels_infanticide)

df_anova$infanticide <- factor(df_anova$infanticide,
                           levels = sort(unique(df_anova$infanticide)),
                           labels = labels_infanticide)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$infanticide

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
c_infanticide <- aov.phylo(dimor ~ soc, tree_c)
summary(c_infanticide)

c = ggplot(df_anova_c, aes(x = infanticide, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Infanticide",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_INFANTICIDE(0.579).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução


df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$infanticide

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
p_infanticide <- aov.phylo(dimor ~ soc, tree_p)
summary(p_infanticide)

p = ggplot(df_anova_p, aes(x = infanticide, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Infaticide",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p


# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_INFANTICIDE(0.058).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução


# ANOVA IN MULTILEVEL SOCIETY AND SD
df_anova <- data.frame(
  genus = df$mds.genus,
  parvorder = df$mds.dados.PARVORDER,
  society = df$mds.dados.MULTILEVEL_SOCIETY,
  dimorfismo = df$align.normas
)

labels_society = c("ABSENCE", "PRESENCE")

df_anova$society <- factor(df_anova$society,
                               labels = labels_society)

df_anova$society <- factor(df_anova$society,
                               levels = sort(unique(df_anova$society)),
                               labels = labels_society)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$genus), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$society

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
c_multilevel <- aov.phylo(dimor ~ soc, tree_c)
summary(c_multilevel)

c = ggplot(df_anova_c, aes(x = society, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of MultiLevel Society",
       y = "Sexual Dimorphism") +
    theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )
c

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/C_SD_MULTILEVEL(0.371).png", plot = c,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução


df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$genus), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$society

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
p_multilevel <- aov.phylo(dimor ~ soc, tree_p)
summary(p_multilevel)

p = ggplot(df_anova_p, aes(x = society, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Multilevel Society",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )
p

# ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/P_SD_MULTILEVEL(0.67).png", plot = p,
#        width = 12,    # largura em inches
#        height = 10,   # altura em inches
#        dpi = 300)    # resolução

##################################################################################

#-----------------------------------
# 1. Função para extrair resultados
#-----------------------------------
extract_model_info <- function(model, model_name) {
  
  # identificar clado
  clade <- ifelse(grepl("^c_", model_name), "Catarrhini", "Platyrrhini")
  
  # identificar preditor
  predictor <- gsub("^[cp]_", "", model_name)
  predictor <- gsub("_lm_model", "continuous", predictor)
  
  # tabela de ANOVA
  aov_tab <- tryCatch(anova(model), error = function(e) NULL)
  
  if (!is.null(aov_tab)) {
    
    row <- aov_tab[1, ]
    
    data.frame(
      Clade = clade,
      Predictor = predictor,
      F_value = as.numeric(row$`F value`),
      p_value = as.numeric(row$`Pr(>F)`)
    )
    
  } else {
    return(NULL)
  }
}

#-----------------------------------
# 2. Lista de modelos
#-----------------------------------
models <- list(
  c_social, c_matings, c_aggression, c_dominance,
  c_all_male, c_furtive, c_infanticide, c_multilevel, c_lm_model,
  p_social, p_matings, p_aggression, p_dominance,
  p_all_male, p_furtive, p_infanticide, p_multilevel, p_lm_model
)

model_names <- c(
  "c_social", "c_matings", "c_aggression", "c_dominance",
  "c_all_male", "c_furtive", "c_infanticide", "c_multilevel", "c_lm_model",
  "p_social", "p_matings", "p_aggression", "p_dominance",
  "p_all_male", "p_furtive", "p_infanticide", "p_multilevel", "p_lm_model"
)

#-----------------------------------
# 3. Criar tabela consolidada
#-----------------------------------
results_table <- map2_df(models, model_names, extract_model_info) %>%
  
  mutate(
    Predictor = recode(Predictor,
                       "social" = "Social organization",
                       "matings" = "Mating system",
                       "aggression" = "Aggression",
                       "dominance" = "Dominance",
                       "all_male" = "All-male groups",
                       "furtive" = "Furtive copulation",
                       "infanticide" = "Infanticide",
                       "multilevel" = "Multilevel society",
                       "lm_model" = "Prop. M/F"
    )
  ) %>%
  
  mutate(
    Clade = case_when(
      Clade == "Catarrhini" ~ "Catarrhini (N = 27)",
      Clade == "Platyrrhini" ~ "Platyrrhini (N = 16)"
    )
  ) %>%
  
  arrange(Clade, p_value)

#-----------------------------------
# 4. Criar tabela gt
#-----------------------------------
gt_table <- results_table %>%
  
  gt(groupname_col = "Clade") %>%
  
  cols_label(
    Predictor = "Predictor",
    F_value = "F",
    p_value = "p"
  ) %>%
  
  fmt_number(
    columns = c(F_value),
    decimals = 2
  ) %>%
  
  fmt_number(
    columns = c(p_value),
    decimals = 3
  ) %>%
  
  # centralizar tudo
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  
  # negrito + centralizado nos clados
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_text(align = "center")
    ),
    locations = cells_row_groups()
  ) %>%
  
  # itálico para p < 0.05
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = p_value,
      rows = p_value < 0.05
    )
  ) %>%
  
  # negrito para p < 0.01
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_value,
      rows = p_value < 0.01
    )
  ) %>%
  
  tab_header(
    title = md("**Phylogenetic ANOVA results**")
  ) %>%
  
  tab_source_note(
    source_note = md(
      "*F-values and p-values derived from phylogenetic ANOVAs. Italics indicate p < 0.05; bold indicates p < 0.01.*"
    )
  )

gt_table
#-----------------------------------
# 5. Exportar
#-----------------------------------
gtsave(gt_table, 
       "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/phylo_anova_table.png",
       vwidth = 1200,
       vheight = 800,
       zoom = 2)
