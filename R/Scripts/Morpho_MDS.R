# set work directory
setwd("~/Dropbox/Doc/Data")

# install and load packages
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')

# reading and cleaning morphological data table
data.morpho = read_xls("CatarrhiniMorphoData.xls")
data.morpho = data.morpho[!is.na(data.morpho$SPECIES), ] # remove missing data rows
data.morpho = data.morpho[-which(data.morpho$SEX == "sexo"), ]
data.morpho = data.morpho[-which(data.morpho$SEX == "0"), ]

# select and cut the data of interest
cut.data.morpho.1 = data.morpho[,c("FAMILY", "GENUS", "SPECIES", "SEX")]
names.species = apply(cut.data.morpho.1[, 2:3], 1, paste, collapse = "_")
cut.data.morpho.1 = cut.data.morpho.1[,-2:-3]
cut.data.morpho.1 = data.frame(cut.data.morpho.1, names.species)
colnames(cut.data.morpho.1)[3] = "SPECIES"


cut.data.morpho.2 = data.morpho[ ,49:87]
cut.data.morpho = data.frame(cut.data.morpho.1, cut.data.morpho.2)

# fix male and female names
cut.data.morpho$SEX[which(cut.data.morpho$SEX == "?male")] = "male"
cut.data.morpho$SEX[which(cut.data.morpho$SEX == "?female")] = "female"

# MDS for morpho data
d = dist(cut.data.morpho[,4:42])
fit = cmdscale(d, eig = TRUE, k = 2)

DIM_1 = fit$points[,1]
DIM_2 = fit$points[,2]

# Data frame for plot
cut.data.morpho = data.frame(cut.data.morpho, DIM_1, DIM_2)

# Plot
mds_all = ggplot(cut.data.morpho, aes(x = DIM_1, y = DIM_2)) +
  geom_point(aes(colour = DIET, shape = SEX), size = 3.5, alpha = 0.4) +
  stat_ellipse(aes(colour = DIET)) +
  theme_classic() +
  labs(colour = "Diet", shape = "Sex")

mds_all