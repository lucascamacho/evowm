# vector correlation between sexual dimorphism, isometric and module vector
# with PCs. Also get the size of the dimorphism vector
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

# load packages and functions
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(ape)){install.packages("ape"); library(ape)}

geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# load data to correlation
medidas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
ancs = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")

# read all vcv matrices
setwd("~/Dropbox/Doc/Data/genus_vcv/")
temp = list.files(pattern = "*.txt")
vcv = lapply(temp, read.csv, header = FALSE, dec = ".", sep = "\t")
names(vcv)  = gsub(".csv", replacement = "", temp)

wrong <- which(sapply(vcv, ncol) == 1)

for(i in wrong){
  
  file <- temp[i]
  
  # tenta outras combinações
  mat <- read.table(file, header = FALSE, dec = ".")
  
  vcv[[i]] <- as.matrix(mat)
}
vcv <- lapply(vcv, as.matrix)
names(vcv) <- gsub("\\.txt$", "", temp)

matings = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
cat = matings$genus[which(matings$dados.PARVORDER == "Catarrhini")]      
platy = matings$genus[which(matings$dados.PARVORDER == "Platyrrhini")]

#
# read and plot phylo tree
commom <- intersect(names(vcv), matings$genus)
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
genus <- sub("_.*", "", tree$tip.label)
keep <- genus %in% commom
tree2 <- drop.tip(tree, tree$tip.label[!keep])
genus2 <- sub("_.*", "", tree2$tip.label)
keep_one <- !duplicated(genus2)
tree_genus <- drop.tip(tree2, tree2$tip.label[!keep_one])
tree_genus$tip.label <- sub("_.*", "", tree_genus$tip.label)
plot(tree_genus)

# get Ancestors
anc_cata = getMRCA(tree_genus, tree_genus$tip.label[1:21])
anc_plat = getMRCA(tree_genus, tree_genus$tip.label[22:37])

# 
dimor_rel = vector()
dimor_abs = vector()
normas = vector()

align_1 = vector()
align_2 = vector()
align_3 = vector()
align_4 = vector()
align_5 = vector()
align_6 = vector()
align_7 = vector()
align_8 = vector()
for(i in 1:length(matings$genus)){
  # chose species
  sp = matings$genus[[i]]
  med = medidas$ByTrait_Averages[[sp]]
  covar = vcv[[sp]]
  
  # se a especie for catarrrhini
  if(sp %in% cat){
    dimor_rel[i] = geomean(med$Fêmeas) / geomean(med$Machos)
    dimorf = unlist(med$Machos) - unlist(med$Fêmeas)
    tam_cra = geomean(med$Machos)
    size = (geomean(med$Machos) + geomean(med$Fêmeas)) / 2
    dimor_abs[i] = mean(dimorf) / size
    
    pc_1 = ancs$PCs[[as.character(anc_cata)]]$PC1
    pc_2 = ancs$PCs[[as.character(anc_cata)]]$PC2
    pc_3 = ancs$PCs[[as.character(anc_cata)]]$PC3
    pc_4 = ancs$PCs[[as.character(anc_cata)]]$PC4
    pc_5 = ancs$PCs[[as.character(anc_cata)]]$PC5
    pc_6 = ancs$PCs[[as.character(anc_cata)]]$PC6
    pc_7 = ancs$PCs[[as.character(anc_cata)]]$PC7
    pc_8 = ancs$PCs[[as.character(anc_cata)]]$PC8
    
    # norm and alignment
    # Male bigger = norm is negative, Female is bigger = norm is positive
    normas_dimorf = norma(dimorf / size)
    #normas[i] = norma(dimorf / size)
    normas[i] = abs(normas_dimorf) * sign(dimor_rel[i] - 1)
    #normas[i] = dimor_rel[i]
    
    align_1[i] = abs(corVector(unlist(pc_1), dimorf / size))
    align_2[i] = abs(corVector(unlist(pc_2), dimorf / size))
    align_3[i] = abs(corVector(unlist(pc_3), dimorf / size))
    align_4[i] = abs(corVector(unlist(pc_4), dimorf / size))
    align_5[i] = abs(corVector(unlist(pc_5), dimorf / size))
    align_6[i] = abs(corVector(unlist(pc_6), dimorf / size))
    align_7[i] = abs(corVector(unlist(pc_7), dimorf / size))
    align_8[i] = abs(corVector(unlist(pc_8), dimorf / size))
    
  }
  
  # e se for platyrrhini
  else{
    dimor_rel[i] = geomean(med$Fêmeas) / geomean(med$Machos)
    dimorf = unlist(med$Machos) - unlist(med$Fêmeas)
    tam_cra = geomean(med$Machos)
    size = (geomean(med$Machos) + geomean(med$Fêmeas)) / 2
    dimor_abs[i] = mean(dimorf) / size

    pc_1 = ancs$PCs[[as.character(anc_plat)]]$PC1
    pc_2 = ancs$PCs[[as.character(anc_plat)]]$PC2
    pc_3 = ancs$PCs[[as.character(anc_plat)]]$PC3
    pc_4 = ancs$PCs[[as.character(anc_plat)]]$PC4
    pc_5 = ancs$PCs[[as.character(anc_plat)]]$PC5
    pc_6 = ancs$PCs[[as.character(anc_plat)]]$PC6
    pc_7 = ancs$PCs[[as.character(anc_plat)]]$PC7
    pc_8 = ancs$PCs[[as.character(anc_plat)]]$PC8
    
    # norm and alignment
    # Male bigger = norm is negative, Female is bigger = norm is positive
    normas_dimorf = norma(dimorf / size)
    #normas[i] = norma(dimorf / size)
    normas[i] = abs(normas_dimorf) * sign(dimor_rel[i] - 1)
    #normas[i] = dimor_rel[i]
    
    align_1[i] = abs(corVector(unlist(pc_1), dimorf / size))
    align_2[i] = abs(corVector(unlist(pc_2), dimorf / size))
    align_3[i] = abs(corVector(unlist(pc_3), dimorf / size))
    align_4[i] = abs(corVector(unlist(pc_4), dimorf / size))
    align_5[i] = abs(corVector(unlist(pc_5), dimorf / size))
    align_6[i] = abs(corVector(unlist(pc_6), dimorf / size))
    align_7[i] = abs(corVector(unlist(pc_7), dimorf / size))
    align_8[i] = abs(corVector(unlist(pc_8), dimorf / size))
  }
}

atuais = data.frame(matings$genus, dimor_abs, normas, align_1, align_2, align_3, align_4, align_5, align_6, align_7, align_8, matings$dados.PARVORDER)

# save
write.csv(atuais, "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv")
