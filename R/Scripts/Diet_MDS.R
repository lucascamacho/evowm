# setwd and load packages
setwd("~/Dropbox/Doc/Data/diets")

# install and load packages
if (!require('MASS')) install.packages('MASS'); library('MASS')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('magrittr')) install.packages('magrittr'); library('magrittr')
if (!require('ade4')) install.packages('ade4'); library('ade4')
if (!require('rgl')) install.packages('rgl'); library('rgl')

# reading diet table
data.diet = read.csv("TrophicGuild.csv")

# species with missing data
#list_species_bad_data = c("Allocebus_trichotis", "Callicebus_nigrifrons", "Cebus_cuscinus", "Cercopithecus_albogularis",
#                          "Cercopithecus_lowei", "Chiropotes_sagulatus", "Daubentonia_madagascariensis", "Eulemur_fulvus", 
#                          "Indri_indri", "Lemur_catta", "Lophocebus_ugandae", "Piliocolobus_rufomitratus",
#                          "Piliocolobus_tephrosceles", "Propithecus_coquereli", "Propithecus_tattersalli", 
#                          "Semnopithecus_schistaceus", "Semnopithecus_vetulus","Semnopithecus_vetulus")


data.diet$TrophicGuild[which(data.diet$TrophicGuild == "Folivore-frugivore")] = 
  "Folivore_frugivore" # fix typos
data.diet = data.diet[!is.na(data.diet$TrophicGuild), ] # remove sp with no class
data.diet = data.diet[!data.diet$TrophicGuild == "", ]

# Diet-similarity matrix
diet = data.diet[,8:15]
diet = lapply(diet, as.numeric)

diet = data.frame(Fruit = diet$Fruit, Leaves = diet$Leaves, Flowers = diet$Flowers, 
                  Seeds = diet$Seeds, Animal = diet$Animal, Gums = diet$Gums, 
                  Nectar = diet$Nectar, Other = diet$Other)

diet = diet / 100 # % in frequency
index = is.na(diet)
diet[index] = 0

rowSums(diet) # not every line sums 1

# normalize diet matrix
diet_norm = t(apply(diet, 1, function(x) x/sum(x)))

# get final data.frame and remove NaN
diet = cbind(data.diet[,1:7], diet_norm)
diet = na.omit(diet)

# all possible combinations between rows
combinations = combn(x = seq_len(length.out = nrow(x = diet)),
                     m = 2, simplify = FALSE)

# create empty matrix similarity
matrix.similarity = matrix(NA, nrow = nrow(diet), ncol = nrow(diet))
colnames(matrix.similarity) = diet$Species
rownames(matrix.similarity) = diet$Species
diag(matrix.similarity) = 1

for(i in 1:length(combinations)){
  #define rows
  row.1 = diet[combinations[[i]][1], ]
  row.2 = diet[combinations[[i]][2], ]
  
  # Calculate diet similarity
  similarity = sqrt(row.1$Fruit * row.2$Fruit) + sqrt(row.1$Leaves * row.2$Leaves) +
    sqrt(row.1$Flowers * row.2$Flowers) + sqrt(row.1$Seeds * row.2$Seeds) +
    sqrt(row.1$Animal * row.2$Animal) + sqrt(row.1$Gums * row.2$Gums) +
    sqrt(row.1$Nectar * row.2$Nectar) + sqrt(row.1$Other * row.2$Other)
  
  # insert similarity in matrix
  matrix.similarity[combinations[[i]][1], combinations[[i]][2]] = similarity
}

matrix.similarity[lower.tri(matrix.similarity)] = t(matrix.similarity)[lower.tri(matrix.similarity)]

# MDS in Diet-similarity matrix
fit = cmdscale(matrix.similarity, eig = TRUE, k = 2)
x = fit$points[,1]
y = fit$points[,2]

# Plot the MDS 
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2",
     main = "Nonmetric MDS", type = "n")
text(x, y, labels = row.names(matrix.similarity), cex = 1)