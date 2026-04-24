setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(readxl)
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

#a_b = d[which(d$nome_especie == 'Alouatta_belzebul'),]
#a_b_m = a_b[which(a_b$Observações == 'Male'),]
#summary(a_b_m)

#a_b = d[which(d$nome_especie == 'Alouatta_belzebul'),]
#a_b_f = a_b[which(a_b$Observações == 'Female'),]
#summary(a_b_f)

#a_f = d[which(d$nome_especie == 'Alouatta_fusca'),]
#a_f_m = a_f[which(a_f$Observações == 'Male'),]
#summary(a_f_m)

#a_f = d[which(d$nome_especie == 'Alouatta_fusca'),]
#a_f_f = a_f[which(a_f$Observações == 'Female'),]
#summary(a_f_f)

#c_j = d[which(d$nome_especie == 'Callithrix_jacchus'),]
#c_j_m = c_j[which(c_j$Observações == 'Male'),]
#summary(c_j_m)

#c_j = d[which(d$nome_especie == 'Callithrix_jacchus'),]
#c_j_f = c_j[which(c_j$Observações == 'Female'),]
#summary(c_j_f)

#c_p = d[which(d$nome_especie == 'Callithrix_penicillata'),]
#c_p_m = c_p[which(c_p$Observações == 'Male'),]
#summary(c_p_m)

#c_p = d[which(d$nome_especie == 'Callithrix_penicillata'),]
#c_p_f = c_p[which(c_p$Observações == 'Female'),]
#summary(c_p_f)

##
library(dplyr)

tabela_medias <- d %>%
  dplyr::group_by(nome_especie, Bioma) %>%
  dplyr::summarise(
    mean_mg       = mean(mg, na.rm = TRUE),
    mean_perc_C   = mean(perc_C, na.rm = TRUE),
    mean_perc_N   = mean(perc_N, na.rm = TRUE),
    mean_Prop_C_N = mean(Prop_C_N, na.rm = TRUE),
    n             = dplyr::n(),
    .groups = "drop"
  )

tabela_medias

tabela_medias <- tabela_medias %>%
  mutate(
    across(starts_with("mean"), ~ round(.x, 2))
  )


tabela_medias <- tabela_medias %>%
  mutate(
    mean_mg       = round(mean_mg, 3),
    mean_perc_C   = round(mean_perc_C, 2),
    mean_perc_N   = round(mean_perc_N, 2),
    mean_Prop_C_N = round(mean_Prop_C_N, 2)
  )

# Biplot
p = ggplot(data = d) +
  #geom_text(aes(x = d13C, y = d15N, color = Bioma, label = nome_especie)) +
  geom_point(aes(x = d13C, y = d15N, color = Bioma, shape = nome_especie), size = 3) + 
  stat_ellipse(aes(x = d13C, y = d15N, group = Bioma, color = Bioma), level = 0.95, linewidth = 1) +
  labs(x = expression(delta^13*C~('\u2030')),
       y = expression(delta^15*N~('\u2030')),
       color = "Biomes",
       shape = "Species") +
  theme_bw()
p

ggsave(filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/Biplot_N_C.png", p,
       dpi = 600, width = 25, height = 15, units = "cm",  bg = "transparent")

# grid de plots 1 por especie e separação M e F
c = ggplot(d, aes(x = d13C, y = d15N, shape = Observações, color = `Tipo da amostra (órgão)`)) +
  geom_point(size = 3, alpha = 0.8) +
  #geom_text(label = d$State) +
  # elipses separando machos/fêmeas dentro de cada espécie
  #stat_ellipse(aes(group = Bioma, color = Bioma), level = 0.95, linewidth = 1) +
  # cria grid 2x2 por espécie
  facet_wrap(~ nome_especie, ncol = 2) +
  labs(
    x = expression(delta^13*C~('\u2030')),
    y = expression(delta^15*N~('\u2030')),
    shape = "Sex",
    color = "Biomes"
  ) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18),     # nome do eixo X
    axis.title.y = element_text(size = 18),     # nome do eixo Y
    strip.text = element_text(face = "bold", size = 16),
    legend.title = element_text(size = 16),     # título da legenda
    legend.text  = element_text(size = 14))      # texto da legenda)
c

ggsave(filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/Biplot_Species.png", c, 
       dpi = 600, width = 25, height = 15, units = "cm",  bg = "transparent")

# So skeleton
d_novo = d[which(d$`Tipo da amostra (órgão)` == 'skeleton'),]




# CALLITHRIX JACCHUS
# Define layout: 1 linha, 2 colunas
#png("~/Dropbox/Doc/Code/evowm/R/Outputs/DensityPlot_Isotopes_Diet_Callithrix.png", width=1600, height=800)
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
       col=c("red","purple"), lwd=2, cex=0.8, bty = 'n')

# --- δ13C ---
d_2 = d[81:120,]
dens_sp_1 <- density(d_2[which(d_2$State == "Ceará"),]$d13C)
dens_sp_2 <- density(d_2[which(d_2$State != "Ceará"),]$d13C)

y_max <- max(dens_sp_1$y, dens_sp_2$y)

plot(dens_sp_1, type="l", xlab="Value", ylab="Density",
     col="red", ylim=c(0, y_max*1.1),
     main="Distribuição de δ13C para os dois grupos")
lines(dens_sp_2, col="purple", lwd=2)
legend("topright", 
       legend=c("C. jacchus Ceará","C. jacchus M.A."),
       col=c("red","purple"), lwd=2, cex=0.8, bty = 'n')
#dev.off()

# ALOUATTA
dens_sp_1 = density(d[1:40,]$d15N)
dens_sp_2 = density(d[41:80,]$d15N)

y_max <- max(dens_sp_1$y, dens_sp_2$y)

#png("~/Dropbox/Doc/Code/evowm/R/Outputs/DensityPlot_Isotopes_Diet_Alouatta.png", width=1600, height=800)
par(mfrow = c(1, 2))  # ou c(2,1) para um em cima do outro
plot(dens_sp_1, type="l", xlab="Value", ylab="Count estimate", col = "red", ylim = c(0, y_max * 1.1),
     main = "Distribuição de δ15N para os dois grupos",
)
lines(dens_sp_2,
      col  = "purple",
      lwd = 2)

# opcional: legenda
legend(
  "topright",
  legend = c(expression(italic("A. belzebul")),
             expression(italic("A. fusca"))),
  col    = c("red", "purple"),
  lwd    = 1,
  cex    = 1,
  bty    = "n")   # remove a caixa

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
legend(
  "topright",
  legend = c(expression(italic("A. belzebul")),
             expression(italic("A. fusca"))),
  col    = c("red", "purple"),
  lwd    = 1,
  cex    = 1,
  bty    = "n")   # remove a caixa
#dev.off()
