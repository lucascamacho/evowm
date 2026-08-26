setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ape)
library(ggplot2)
library(tidyverse)
library(ggridges)

load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/vcv.RData")

load("~/Desktop/Primaset_RawData.RData")

load("~/Desktop/Primatrees.RData")
species <- names(vcv)
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

medias <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData")
sd <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData")

clades <- read.csv("~/Desktop/taxonomy.csv")

#
measure_info <- data.frame(
  Measure = distn,
  Module = c(
    "Face",       # IS-PM
    "Face",       # IS-NSL
    "Face",       # IS-PNS
    "Face",       # PM-ZS
    "Face",       # PM-ZI
    "Face",       # PM-MT
    "Face",       # NSL-NA
    "Face",       # NSL-ZS
    "Face",       # NSL-ZI
    "Neurocranium", # NA-BR
    "Neurocranium", # NA-FM
    "Face",        # NA-PNS
    "Neurocranium",# BR-PT
    "Neurocranium",# BR-APET
    "Neurocranium",# PT-FM
    "Neurocranium",# PT-APET
    "Neurocranium",# PT-BA
    "Neurocranium",# PT-EAM
    "Face",        # PT-ZYGO
    "Both", # PT-TSP
    "Neurocranium",# FM-ZS
    "Neurocranium",# FM-MT
    "Face",        # ZS-ZI
    "Face",        # ZI-MT
    "Face",        # ZI-ZYGO
    "Face",        # ZI-TSP
    "Face",        # MT-PNS
    "Neurocranium",# PNS-APET
    "Neurocranium",# APET-BA
    "Neurocranium",# APET-TS
    "Neurocranium",# BA-EAM
    "Face",        # EAM-ZYGO
    "Face",        # ZYGO-TSP
    "Neurocranium",# LD-AS
    "Neurocranium",# BR-LD
    "Neurocranium",# OPI-LD
    "Neurocranium",# PT-AS
    "Neurocranium",# JP-AS
    "Neurocranium" # BA-OPI
  )
)

sexual_dimorphism_long <- sexual_dimorphism_df %>%
  select(Genus, all_of(as.character(1:39))) %>%
  pivot_longer(
    cols = -Genus,
    names_to = "Measure_ID",
    values_to = "Sexual_Dimorphism"
  ) %>%
  mutate(
    Measure_ID = as.integer(Measure_ID)
  ) %>%
  left_join(
    measure_info %>%
      mutate(Measure_ID = 1:39),
    by = "Measure_ID"
  )

p_sexual_dimorphism <- ggplot(
  sexual_dimorphism_long,
  aes(
    x = Sexual_Dimorphism,
    y = Measure
  )
) +
  
  geom_density_ridges(
    aes(
      height = after_stat(density),
      group = Measure
    ),
    stat = "density",
    scale = 1.2,
    rel_min_height = 0.01,
    alpha = 0.6
  ) +
  
  geom_point(
    position = position_jitter(
      height = 0.08,
      width = 0
    ),
    size = 1.5,
    alpha = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  facet_grid(
    Module ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  labs(
    x = "Log Sexual dimorphism",
    y = NULL
  ) +
  
  theme_classic()

p_sexual_dimorphism

sexual_dimorphism_long <- sexual_dimorphism_long %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Genus" = "GENUS")
  )

p_sexual_dimorphism <- ggplot(
  sexual_dimorphism_long,
  aes(
    x = Sexual_Dimorphism,
    y = Measure
  )
) +
  
  geom_density_ridges(
    aes(
      height = after_stat(density),
      group = Measure
    ),
    stat = "density",
    scale = 1.2,
    rel_min_height = 0.01,
    alpha = 0.5
  ) +
  
  geom_point(
    aes(color = CLADE),
    position = position_jitter(
      height = 0.08,
      width = 0
    ),
    size = 1.7,
    alpha = 0.7
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  facet_grid(
    Module ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  labs(
    x = "Log Sexual dimorphism: (log(M) - log(F))",
    y = NULL,
    color = "Clade"
  ) +
  
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

p_sexual_dimorphism

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Figure_RidgePlot.png",
  plot = p_sexual_dimorphism,
  width = 12,
  height = 6,
  dpi = 300
)
