# set WD and functions
setwd("~/Dropbox/Doc/Output/")

if(!require(ape)){install.packages("ape"); library(ape)}

# load mating systems data 
matings = read.table("~/Dropbox/Doc/Data/wos_mating_systems/primates_mating_systems.csv",
                     sep = ",", header = TRUE)

# get "interested" species
species_tree = vector()
for(i in 1:nrow(matings)){
  if(is.na(matings$SUBSPECIES[i])){
    species_tree[i] = paste(matings$GENUS[i], matings$SPECIES[i], sep="_")}
  else{
    species_tree[i] = paste(matings$GENUS[i], matings$SPECIES[i], matings$SUBSPECIES[i], sep="_")
  }
}

# Read MDS and Dimorps Cor PCS data
mds = read.table("Dimorp_MDS_OrgVars.txt", header = TRUE, row.names = NULL) #AQUI
dimor = read.table("Pcs_Dimorp_Norma.txt", header = TRUE, row.names = NULL, sep = "\t")
mds$species = paste(mds$row.names, mds$species, sep = "_")

index = match(dimor$names.vcv., mds$species)
index = index[!is.na(index)]
mds_new = mds[index,]

index = match(mds_new$species, dimor$names.vcv.)
dimor_new = dimor[index,]

dimor = dimor_new
mds = mds_new
species_tree = dimor$names.vcv.

#read tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
cat.tree = keep.tip(tree, species_tree)

#
co = data.frame()
um = cor(mds$mds1, dimor$dimorp_PC1)
dois = cor(mds$mds1, dimor$dimorp_PC2)
tres = cor(mds$mds1, dimor$dimorp_PC3)
quatro = cor(mds$mds1, dimor$dimorp_PC4)

cinco = cor(mds$mds2, dimor$dimorp_PC1)
seis = cor(mds$mds2, dimor$dimorp_PC2)
sete = cor(mds$mds2, dimor$dimorp_PC3)
oito = cor(mds$mds2, dimor$dimorp_PC4)

co = rbind(co, c(um, dois, tres, quatro, cinco, seis, sete, oito))

for(i in 1:8){
  co[2,i] =  0.5 * ((log(1 + co[1,i]) / log(1 - co[1,i])))
}

pes = vector()
for(i in 1:4){
  mds_1 = setNames(mds$mds1, mds$species)
  dimor_pc1 = setNames(dimor[,i+2], mds$species)
  
  pic_mds = pic(mds_1, cat.tree)
  pic_dimor = pic(dimor_pc1, cat.tree)

  co[3,i] = 0.5 * ((log(1 + cor(pic_mds, pic_dimor)) / log(1 - cor(pic_mds, pic_dimor))))

  z = 0.5 * ((log(1 + cor(pic_mds, pic_dimor)) / log(1 - cor(pic_mds, pic_dimor))))
  zse = 1 / sqrt(length(pic_mds) - 3)
  
  pes[i] = min(pnorm(z, sd=zse), pnorm(z, lower.tail = F, sd = zse))*2
  
}

for(i in 1:4){
  mds_2 = setNames(mds$mds2, mds$species)
  dimor_pc1 = setNames(dimor[,i+2], mds$species)
  
  pic_mds = pic(mds_2, cat.tree)
  pic_dimor = pic(dimor_pc1, cat.tree)
  
  co[3, i+4] = 0.5 * ((log(1 + cor(pic_mds, pic_dimor)) / log(1 - cor(pic_mds, pic_dimor))))
  
  z = 0.5 * ((log(1 + cor(pic_mds, pic_dimor)) / log(1 - cor(pic_mds, pic_dimor))))
  zse = 1/sqrt(71 - 3)
  pes[i+4] = min(pnorm(z, sd=zse), pnorm(z, lower.tail=F, sd=zse))*2
  
}

colnames(co) = c("mds1_pc1", "mds1_pc2", "mds1_pc3", "mds1_pc4", "mds2_pc1", "mds2_pc2", "mds2_pc3", "mds2_pc4")
row.names(co) = c("simple", "Fisher", "PIC")


co = rbind(co, pes)
row.names(co) = c("simple", "Fisher", "PIC", "p-value")

write.table(co, file = "~/Dropbox/Doc/Output/table_corr_pic_z.txt")

### cor MDS Org Var
cor(mds$mds2, mds$SOCIAL_ORGANIZATION)
cor(mds$mds2, mds$MATING_SYSTEM)
cor(mds$mds2, mds$PROP_MALES_FEMALES)
cor(mds$mds2, mds$DOMINANCE)
cor(mds$mds2, mds$AGGRESSION)
