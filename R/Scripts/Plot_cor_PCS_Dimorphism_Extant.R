setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(ggplot2)
library(viridis)
library(dplyr)
library(scales)
library(gridExtra)
library(patchwork)
library(ggpmisc)

#
data = read.csv2("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv", sep = ",")

#
cols <- 3:12

data[ , cols] <- lapply(data[ , cols], function(x) {
  if(!is.numeric(x)) as.numeric(as.character(x)) else x
})

# Criar coluna de grupos
data$grupo <- ifelse(1:nrow(data) <= 38, "Catarrhini", "Platyrrhini")

mycols <- c("#1b9e77", "#7570b3")

# p1 com R²
p1 <- ggplot(data, aes(x = align_1, y = normas, color = grupo, fill = grupo)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.4,
              linetype = "dashed", alpha = 0.2) +
  geom_point(size = 2.5, alpha = 0.9, stroke = 2) +
  stat_poly_eq(
    aes(label = paste(..rr.label..), color = grupo),
    formula = y ~ x,
    parse = TRUE,
    label.x.npc = "left",
    label.y.npc = c(0.95, 0.85),
    size = 6
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  scale_color_manual(values = mycols) +
  scale_fill_manual(values = mycols) +
  theme_classic(base_size = 14) +
  theme_bw(base_size = 14) +
  labs(
    x = "PC1-Sexual Dimorphism Alignments",
    y = "Normalized Magnitude of Sexual Dimorphism",
    color = "Parvorder",
    fill  = "Parvorder"
  ) +
  theme(
    #    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold")
  )

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_Dimorphism_Extant.png", plot = p1,
       width = 12,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


##########################
# All values of aligns with SD and PCs 1 to 8
# Seleciona as colunas
cols <- 5:12

# Define grid 4x2 (4 linhas, 2 colunas)
png(
  filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/Hist_All_cor_PCS_Dimorphism_Extant.png",
  width = 4000,      # largura em pixels
  height = 3000,     # altura em pixels
  res = 300           # resolução em dpi
)
par(mfrow = c(4, 2), mar = c(3, 3, 2, 1))
for(i in cols){
  hist(data[[i]],
       main = colnames(data)[i],
       xlab = "",
       ylab = "",
       col = "lightblue",
       border = "white")
}
dev.off()

# Extract R2
# ajustar os modelos
mods <- data %>%
  group_split(grupo) %>%
  setNames(levels(data$grupo)) %>%
  lapply(function(df) lm(normas ~ align_1, data = df))

# extrair R²
r2_vals <- sapply(mods, function(m) summary(m)$r.squared)
r2_vals

