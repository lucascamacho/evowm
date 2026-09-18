msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")

# check names and remove doubts
msrs$SEX[which(msrs$SEX == "?female")] = "female"
msrs$SEX[which(msrs$SEX == "?male")] = "male"

# remove uncertain sex
msrs$SEX[msrs$SEX == "0"] = NA
msrs$SEX[msrs$SEX == ""] = NA
msrs$SEX[msrs$SEX == "sexo"] = NA

msrs = msrs[complete.cases(msrs$SEX), ]

# new column with species names
species = vector()
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUBSPECIES[i])){
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS[i], msrs$SPECIES[i], msrs$SUBSPECIES[i], sep="_")
  }
}

msrs = cbind(species, msrs)
species = unique(species)

data_museus = data.frame()
for(i in 1:length(species)){
  species_subset = msrs[which(msrs$species == species[i]),]
  museus = unique(species_subset$MUSEUM)
  for(j in 1:length(museus)){
    museu_subset = species_subset[species_subset$MUSEUM == museus[j],]
    
    n_males = length(which(museu_subset$SEX == "male"))
    n_females = length(which(museu_subset$SEX == "female"))
    
    data_final = c(species[i], museus[j], n_males, n_females)
    data_museus = rbind(data_museus, data_final)
  }
}

colnames(data_museus) = c("Species", "Museu", "N_Males", "N_Females")

write.xlsx(data_museus, file = "~/Dropbox/Doc/Code/evowm/R/Novo_Output/Numero_individuos_museus_Catarrhini.xlsx")
