# set work directory
setwd("~/Dropbox/Doc/Data/wos_mating_systems/Catarrhini")

# load packages
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

# read all the species files
temp = list.files(pattern="*.xls")
excells = lapply(temp, read_xls)
names(excells)  = gsub(".xls", replacement= "", temp)

# genus names
gen = unique(sapply(strsplit(temp, "_"), `[`, 1))

# define empty vectors
total = vector() # total number of articles from WoS
choice = vector() # picked articles for mating system research

for(i in 1:length(temp)){
  total[i] = length(pull(excells[[i]][1])) # row number for each species file
  choice[i] = sum(excells[[i]][,1]) # sum the first column for each species file
}

# basic statistics
sum(total)
sum(choice)

# basic plots showing the total and picked number distribution per species
hist(total, main = "Histogram of articles from WoS by species",
     ylab = "Frequency",
     xlab = "Species")
plot(total, main = "Total number of articles from WoS by species",
     ylab = "Number of articles",
     xlab = "Species",
     labels = names(excells))

hist(choice, main = "Histogram of articles picked for Mating System research",
     ylab = "Frequency",
     xlab = "Species")
plot(choice, main = "Total number of articles picked for Mating System research",
     ylab = "Number of articles",
     xlab = "Species")

# Number of articles separating by genus
total_genus = vector()
choice_genus = vector()

for(i in 1:length(gen)){
  temp = list.files(pattern = gen[i]) # find and read the files for each genus
  excells = lapply(temp, read_xls)
  names(excells)  = gsub(".xls", replacement= "", temp)
  
  for(j in 1:length(excells)){
    total_genus[i] = length(pull(excells[[j]][1])) # number of articles of a genus
    choice_genus[i] = sum(as.numeric(pull(excells[[j]][1]))) # picked articles of a genus
  }
}

# Basic statistics and plots separated by genus of total and picked number of articles
hist(total_genus)
plot(total_genus)

hist(choice_genus)
plot(choice_genus)