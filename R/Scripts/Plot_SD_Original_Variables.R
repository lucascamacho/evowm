setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(ape)
library(ggplot2)
library(dplyr)
library(ggpmisc)
library(geiger)

# datas
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

#Load and plot the phylogeny
# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = read.nexus(filename)
species = mds$especies
tree$tip.label[tree$tip.label == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"
tree = drop.tip(tree, setdiff(tree$tip.label, species))
tree_c = drop.tip(tree, setdiff(tree$tip.label, align$matings.especies[1:38]))
tree_p = drop.tip(tree, setdiff(tree$tip.label, align$matings.especies[39:64]))

df = data.frame(mds$especies, mds$dados.PARVORDER, mds$dados.SOCIAL_ORGANIZATION, mds$dados.MATING_SYSTEM,
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
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6, col = "red") +
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

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/F_M_cor_Dimorphism_Extant_integrationIndex.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução

# ANOVA in SOCIAL ORG AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
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
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$especie), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$social
names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

anova_model = aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

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

ggsave("~/Dropbox/Doc/C_SD_SOCIAL_ORG.png", plot = c,
       width = 12,    # largura em inches
       height = 10,   # altura em inches
       dpi = 300)    # resolução

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$social

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

anova_model = aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)

p = ggplot(df_anova_p, aes(x = social, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = " ",
       y = "Magnitude of Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

ggsave("~/Dropbox/Doc/P_SD_SOCIAL_ORG.png", plot = p,
       width = 12,    # largura em inches
       height = 10,   # altura em inches
       dpi = 300)    # resolução

# ANOVA in MATING SYSTEMS AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
  parvorder = df$mds.dados.PARVORDER,
  matings = df$mds.dados.MATING_SYSTEM,
  dimorfismo = df$align.normas
)

labels_matings = c(1,2,3,4)

#labels_matings = c("MONOGAMY", "POLYGYNANDROUS or  MONOGAMY/POLYGYNY or MONOGAMY/POLYANDROUS", "POLYGYNY/POLYGYNANDROUS",
#                   "POLYGYNY or POLYANDROUS")

df_anova$matings <- factor(df_anova$matings,
                          labels = labels_matings)

df_anova$matings <- factor(df_anova$matings,
                          levels = sort(unique(df_anova$matings)),
                          labels = labels_matings)

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$especie), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$matings

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

anova_model = aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

c = ggplot(df_anova_c, aes(x = matings, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Mating Systems",
       y = " ") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c


ggsave("~/Dropbox/Doc/C_SD_MATING_SYSTEM.png", plot = c,
       width = 12,    # largura em inches
       height = 10,   # altura em inches
       dpi = 300)    # resolução

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$matings

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

anova_model = aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)

p = ggplot(df_anova_p, aes(x = matings, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Mating Systems",
       y = " ") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

ggsave("~/Dropbox/Doc/P_SD_MATING_SYSTEM.png", plot = p,
       width = 12,    # largura em inches
       height = 10,   # altura em inches
       dpi = 300)    # resolução

# LM in PROP M/F AND SD
df_lm <- data.frame(
  especie = df$mds.especies,
  parvorder = df$mds.dados.PARVORDER,
  prop = df$mds.dados.PROP_MALES_FEMALES,
  dimorfismo = df$align.normas
)

df_lm_c = df_lm[which(df_lm$parvorder == "Catarrhini"),]
lm_model <- lm(dimorfismo ~ prop, data = df_lm_c)
r2 <- summary(lm_model)$r.squared %>% round(3)
summary(lm_model)

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

df_lm_p = df_lm[which(df_lm$parvorder == "Platyrrhini"),]
lm_model <- lm(dimorfismo ~ prop, data = df_lm_p)
r2 <- summary(lm_model)$r.squared %>% round(3)
summary(lm_model)

p = ggplot(df_lm_p, aes(x = prop, y = dimorfismo)) +
  geom_point(size = 4, alpha = 0.9, stroke = 2, color = "#7570b3") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
  labs(
    x = "Proportion M/F",
    y = " "
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

ggsave("~/Dropbox/Doc/P_SD_PROP.png", plot = p,
       width = 12,    # largura em inches
       height = 10,   # altura em inches
       dpi = 300)    # resolução


# ANOVA in DOMINANCE AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
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
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$especie), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$dominance

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

c = ggplot(df_anova_c, aes(x = dominance, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Dominance",
       y = " ") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$dominance

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)

p = ggplot(df_anova_p, aes(x = dominance, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Dominance",
       y = " ") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

# ANOVA in AGGRESSION AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
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
df_anova_c = df_anova_c[match(tree_c$tip.label, df_anova_c$especie), ]

dimor = df_anova_c$dimorfismo
soc = df_anova_c$aggression

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

c = ggplot(df_anova_c, aes(x = aggression, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Aggression",
       y = " ") +
  theme_classic() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
c

ggsave("~/Dropbox/Doc/C_SD_AGGRESSION.png", plot = c,
       width = 12,    # largura em inches
       height = 10,   # altura em inches
       dpi = 300)    # resolução


df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$aggression

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)

p = ggplot(df_anova_p, aes(x = aggression, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 4, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = " ",
       y = " ") +
  theme_classic() +
  theme(
    #axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
p

ggsave("~/Dropbox/Doc/P_SD_AGGRESSION.png", plot = p,
       width = 12,    # largura em inches
       height = 10,   # altura em inches
       dpi = 300)    # resolução

# ANOVA in ALL MALE GROUPS AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
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

dimor = df_anova_c$dimorfismo
soc = df_anova_c$all_male

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

c = ggplot(df_anova_c, aes(x = all_male, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
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
c

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$all_male

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)


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

# ANOVA in FURTIVE COPULATION AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
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

dimor = df_anova_c$dimorfismo
soc = df_anova_c$furtive

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

c = ggplot(df_anova_c, aes(x = furtive, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Furtive Copulations",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )
c

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$furtive

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)

p = ggplot(df_anova_p, aes(x = furtive, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Furtive Copulation",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )
p

# ANOVA IN INFANTICIDE AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
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

dimor = df_anova_c$dimorfismo
soc = df_anova_c$infanticide

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

c = ggplot(df_anova_c, aes(x = infanticide, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#1b9e77", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Infanticide",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )
c

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$infanticide

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)

p = ggplot(df_anova_p, aes(x = infanticide, y = dimorfismo)) +
  geom_jitter(width = 0.05, color = "#7570b3", size = 2.5, alpha = 0.9, stroke = 2) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  labs(x = "Evidences of Infaticide",
       y = "Sexual Dimorphism") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )
p


# ANOVA IN MULTILEVEL SOCIETY AND SD
df_anova <- data.frame(
  especie = df$mds.especies,
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

dimor = df_anova_c$dimorfismo
soc = df_anova_c$society

names(dimor) = tree_c$tip.label
names(soc) = tree_c$tip.label

df_anova_c = df_anova[which(df_anova$parvorder == "Catarrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_c)
summary(anova_model)

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

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
df_anova_p = df_anova_p[match(tree_p$tip.label, df_anova_p$especie), ]

dimor = df_anova_p$dimorfismo
soc = df_anova_p$society

names(dimor) = tree_p$tip.label
names(soc) = tree_p$tip.label

df_anova_p = df_anova[which(df_anova$parvorder == "Platyrrhini"),]
anova_model <- aov.phylo(dimor ~ soc, tree_p)
summary(anova_model)

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

