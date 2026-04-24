setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

library(ggplot2)
library(viridis)
library(dplyr)
library(scales)
library(gridExtra)
library(patchwork)
library(ggpmisc)
library(tidyr)

#
data = read.csv2("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv", sep = ",")

#
cols <- 3:12

data[ , cols] <- lapply(data[ , cols], function(x) {
  if(!is.numeric(x)) as.numeric(as.character(x)) else x
})

# Criar coluna de grupos
mycols <- c("#1b9e77", "#7570b3")

# p1 com R²
p1 <- ggplot(data, aes(x = atanh(align_1), y = normas, color = matings.dados.PARVORDER, fill = matings.dados.PARVORDER)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.4,
              linetype = "dashed", alpha = 0.2) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6, col = "red") +
  stat_poly_eq(
    aes(label = paste(..rr.label..), color = matings.dados.PARVORDER),
    formula = y ~ x,
    parse = TRUE,
    size = 6,
    label.x = "left",
    label.y = "bottom"
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  scale_color_manual(values = mycols) +
  scale_fill_manual(values = mycols) +
  theme_classic(base_size = 14) +
  theme_bw(base_size = 14) +
  labs(
    x = "Pmax-Sexual Dimorphism Alignments",
    y = "Sexual Dimorphism",
    color = "Parvorder",
    fill  = "Parvorder"
  ) +
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

p1

#summary(lm(normas ~ atanh(align_1), data))

ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/F_M_cor_PCS_Dimorphism_Extant.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução

##########################
# All values of aligns with SD and PCs 1 to 8
# Seleciona as colunas
cols <- 5:8
titles <- paste("SD and PC", 1:4)

df_long <- data %>%
  dplyr::select(all_of(cols)) %>%
  setNames(titles) %>%
  tidyr::pivot_longer(
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
    fill = "#6baed6",
    color = "white",
    alpha = 0.8
  ) +
  geom_density(
    color = "#08519c",   # azul mais escuro
    linewidth = 1.2,
    alpha = 0.7
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
    title = "Alignment Between Sexual Dimorphism and Principal Components"
  )
p

# Salva o gráfico em alta resolução
ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Hist_All_cor_PCS_Dimorphism_Extant.png",
  plot = p,
  width = 12,
  height = 7,
  dpi = 300
)


# Extract R2 de Platyrrhini e Catarrhini
# ajustar os modelos
mods <- data %>%
  group_split(matings.dados.PARVORDER) %>%
  setNames(levels(data$matings.dados.PARVORDER)) %>%
  lapply(function(df) lm(normas ~ atanh(align_1), data = df))

# extrair R²
r2_vals <- sapply(mods, function(m) summary(m)$r.squared)
r2_vals

