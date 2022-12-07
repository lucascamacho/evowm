# Initial exploratory analisys in OWM data
# setwd and load packages
setwd("~/Dropbox/Doutorado/Data/")
library(readxl)
library(ggplot2)
install.packages("ggpubr")
library(ggpubr)

# read table
data = read_excel("Catarrhini - BAOPI por último - médias.xls")
head(data)

only_measures = data[,49:87]

genus = data.frame(data[,8], data[,49:87])

aggregate(genus[,2:40], list(genero = genus[,1]), mean)

ggplot(melt(genus), aes(value, group = GENUS, fill = GENUS)) + 
  geom_histogram() + 
  facet_wrap(~variable, scale = "free")

ggplot(genus, aes(genus$ISPN, genus$ISNSL, group = GENUS, color = GENUS)) + 
  geom_point() +
  geom_smooth(method = "lm", aes(group = GENUS))


mean(only_measures)


cov(only_measures)

library(ggplot2)
library(reshape2)
if(!require(viridis)) install.packages("viridis")
library(viridis)
plotMatrix <- function (corMat, file = NULL) {
  diag(corMat) <- NA
  n_traits = nrow(corMat)
  myPalette <- viridis(50)
  ## Se quiser uma paleta All American, use essa linha em vez da anterior
  #myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 50)
  m.rs = melt(corMat)
  m.rs$Var1 <- factor(m.rs$Var1, levels = m.rs$Var1[n_traits:1])
  m.rs.position = m.rs
  m.rs.position$Var1 <- as.numeric(m.rs.position$Var1)
  m.rs.position$Var2 <- as.numeric(m.rs.position$Var2)
  m.rs.position$value= round(m.rs.position$value, 2)
  m.rs.position$value[is.na(m.rs.position$value)] <- levels(m.rs$Var1)[n_traits:1]
  p <-
    ggplot (m.rs) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradientn(name = '', colours = myPalette, limits = c(-0.15, 1)) +
    labs(x = NULL, y = NULL) +
    geom_text(data = m.rs.position, aes(x = Var2, y = Var1, label = value)) +
    theme_bw()
  if(!is.null(file)) cowplot::save_plot(plot = p, file)
  return(p)
}

cor_A = cor(only_measures)
plotMatrix(cor_A)

pcs = eigen(cov(only_measures))

plot(pcs$vectors[,2])

hist(data[,49])

ggdensity(only_measures$ISNSL,
          main = "density plot",
          xlab = "Values")
sapply(only_measures[1:4900,], shapiro.test)

medias = sapply(only_measures, mean)
desvio = sapply(only_measures, sd)

mean(medias)
c

hist(only_measures[,1])
str(only_measures)


mean_plot = ggplot(data = only_measures) +
  geom_point(aes())






# basic statistics for each cranial measure
m = sapply(data[,8:46], mean)
dp = sapply(data[,8:46], sd)

# 
per_species = aggregate(data[,8:46], list(sex = data$SEX), mean)


# Função que calcula o coeficiente de variação
cv = function(x) sd(x)/mean(x)

# Coeficiente de variação
cef = sapply(data[,8:46], cv)

