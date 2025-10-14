setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(openxlsx)
library(dplyr)
library(ggplot2)

d = read_xls("~/Dropbox/Doc/Data/diet/local_completo__Res Finais PLA 1_167 Plinio_Lucas_macaco.xls")
d = d[-168, -c(18, 19)]

d <- d %>%
  dplyr::rename(
    d13C = '13C', 
    d15N = '15N',
    perc_C = '%C',
    perc_N = '%N',
    Prop_C_N = 'C/N'
  ) %>%
  mutate(nome_especie = paste(...5, ...6, sep = "_"))

# Biplot
p = ggplot(data = d) +
  geom_text(aes(x = d13C, y = d15N, color = Bioma, label = nome_especie)) +
  stat_ellipse(aes(x = d13C, y = d15N, group = Bioma, color = Bioma), level = 0.95, linewidth = 1) +
  labs(x = expression(delta^13*C~('\u2030')),
       y = expression(delta^15*N~('\u2030'))) +
  theme_bw()

ggsave(p, filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/Biplot_N_C.pdf", dpi = 600,
       width = 35, height = 15, units = "cm",  bg = "transparent")

# grid de plots 1 por especie e separação M e F
c = ggplot(d, aes(x = d13C, y = d15N, shape = Observações, color = Observações)) +
  geom_point(size = 3, alpha = 0.8) +
  #geom_text(label = d$State) +
  # elipses separando machos/fêmeas dentro de cada espécie
  stat_ellipse(aes(group = Observações, color = Observações), level = 0.95, linewidth = 1) +
  # cria grid 2x2 por espécie
  facet_wrap(~ nome_especie, ncol = 2) +
  labs(
    x = expression(delta^13*C~('\u2030')),
    y = expression(delta^15*N~('\u2030')),
    shape = "Sexo",
    color = "Sexo"
  ) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 12))

ggsave(c, filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/Biplot_Species.pdf", dpi = 600,
       width = 65, height = 35, units = "cm",  bg = "transparent")


# CALLITHRIX JACCHUS
# Define layout: 1 linha, 2 colunas
png("~/Dropbox/Doc/Code/evowm/R/Outputs/DensityPlot_Isotopes_Diet_Callithrix.png", width=1600, height=800)
par(mfrow = c(1, 2))  # ou c(2,1) para um em cima do outro

# --- δ15N ---
d_1 = d[81:120,]
dens_sp_1 <- density(d_1[which(d_1$State == "Ceará"),]$d15N)
dens_sp_2 <- density(d_1[which(d_1$State != "Ceará"),]$d15N)

y_max <- max(dens_sp_1$y, dens_sp_2$y)

plot(dens_sp_1, type="l", xlab="Value", ylab="Density",
     col="red", ylim=c(0, y_max*1.1),
     main="Distribuição de δ15N para os dois grupos")
lines(dens_sp_2, col="purple", lwd=2)
legend("topright", legend=c("C. jacchus Ceará","C. jacchus M.A."),
       col=c("red","purple"), lwd=2, cex=0.8)

# --- δ13C ---
d_2 = d[81:120,]
dens_sp_1 <- density(d_2[which(d_2$State == "Ceará"),]$d13C)
dens_sp_2 <- density(d_2[which(d_2$State != "Ceará"),]$d13C)

y_max <- max(dens_sp_1$y, dens_sp_2$y)

plot(dens_sp_1, type="l", xlab="Value", ylab="Density",
     col="red", ylim=c(0, y_max*1.1),
     main="Distribuição de δ13C para os dois grupos")
lines(dens_sp_2, col="purple", lwd=2)
legend("topright", legend=c("C. jacchus Ceará","C. jacchus M.A."),
       col=c("red","purple"), lwd=2, cex=0.8)
dev.off()

# ALOUATTA
dens_sp_1 = density(d[1:40,]$d15N)
dens_sp_2 = density(d[41:80,]$d15N)

y_max <- max(dens_sp_1$y, dens_sp_2$y)

png("~/Dropbox/Doc/Code/evowm/R/Outputs/DensityPlot_Isotopes_Diet_Alouatta.png", width=1600, height=800)
par(mfrow = c(1, 2))  # ou c(2,1) para um em cima do outro
plot(dens_sp_1, type="l", xlab="Value", ylab="Count estimate", col = "red", ylim = c(0, y_max * 1.1),
     main = "Distribuição de δ15N para os dois grupos",
)
lines(dens_sp_2,
      col  = "purple",
      lwd = 2)

# opcional: legenda
legend("topright",
       legend = c("A. belzebul", "A. fusca"),
       col    = c("red", "purple"),
       lwd    = 2)


dens_sp_1 = density(d[1:40,]$d13C)
dens_sp_2 = density(d[41:80,]$d13C)

y_max <- max(dens_sp_1$y, dens_sp_2$y)

plot(dens_sp_1, type="l", xlab="Value", ylab="Count estimate", col = "red", ylim = c(0, y_max * 1.1),
     main = "Distribuição de δ13C para os dois grupos",
)
lines(dens_sp_2,
      col  = "purple",
      lwd = 2)

# opcional: legenda
legend("topright",
       legend = c("A. belzebul", "A. fusca"),
       col    = c("red", "purple"),
       lwd    = 2)
dev.off()
