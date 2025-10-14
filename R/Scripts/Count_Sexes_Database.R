setwd("~/Dropbox/Doc/Code/evowm/R/Scripts")

# Read and unique species
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")

# museums
mus = msrs[which(msrs$MUSEUM == ""),]

# insert NA when there is no SUBSPECIES name
msrs$SUBSPECIES[which(msrs$SUBSPECIES == "0")] = NA
msrs$SUBSPECIES[which(msrs$SUBSPECIES == "")] = NA

# check names and remove doubts
msrs$SEX[which(msrs$SEX == "?female")] = "female"
msrs$SEX[which(msrs$SEX == "?male")] = "male"

# remove incertain sex
msrs$SEX[msrs$SEX == "0"] = NA
msrs$SEX[msrs$SEX == ""] = NA
msrs$SEX[msrs$SEX == "sexo"] = NA

msrs = msrs[complete.cases(msrs$SEX), ]

species = vector()
# new column with species names
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUBSPECIES[i])){
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], msrs$SUBSPECIES[i], sep="_")
  }
}

msrs = cbind(species, msrs)
species = unique(species)

count = data.frame()
for(i in 1:length(species)){
  # choose species
  sub_species = msrs[which(msrs$species == species[i]), ]
  
  males = table(sub_species$SEX)[2]
  females = table(sub_species$SEX)[1]
  
  count = rbind(count, c(males, females))
}  
  
# make data frame and save
count = cbind(species, count)
write.table(count, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/Species_MF_Count_Database.csv")
