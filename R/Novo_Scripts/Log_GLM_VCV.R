# set WD and functions
setwd("~/Dropbox/Doc/Output/")

# load packages
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

# Read and unique species
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv", dec = ",", sep = ",")

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

###############################################################################################################

# Allenopithecus nigroviridis
nigrovir = msrs[which(msrs$SPECIES == "nigroviridis"),]

fit = lm(log(as.matrix(nigrovir[,50:88])) ~ nigrovir$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Allenopithecus_nigroviridis.txt", col.names = NA)

# Bunopithecus hoolok
hoolok = msrs[which(msrs$SPECIES == "hoolock"),]
hoolok = hoolok[-c(51,52),]
index = which(hoolok$SUBSPECIES == "")
hoolok$SUBSPECIES[index] = "hoolock"

fit = lm(log(as.matrix(hoolok[,50:88])) ~ hoolok$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Bunopithecus_hoolock.txt", col.names = NA)

# Cercocebus agilis
agilis = msrs[which(msrs$SPECIES == "agilis"),]
agilis = agilis[which(agilis$GENUS == "Cercocebus"),]

fit = lm(log(as.matrix(agilis[,50:88])) ~ agilis$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercocebus_agilis.txt", col.names = NA)

# Cercocebus atys
atys = msrs[which(msrs$SPECIES == "atys"),]
atys$SUBSPECIES[which(atys$SUBSPECIES == "")] = "atys"
atys$SUBSPECIES[which(atys$SUBSPECIES == "0")] = "atys"

fit = manova(log(as.matrix(atys[,50:88])) ~ atys$SEX + atys$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercocebus_atys.txt", col.names = NA)

# Cercocebus chrysogaster
chryso = msrs[which(msrs$SPECIES == "chrysogaster"),]
# grupo irmao
agi = msrs[which(msrs$SPECIES == "agilis"),]
agi = agi[which(agi$GENUS == "Cercocebus"),]
gale = msrs[which(msrs$SPECIES == "galeritus"),]

chryso = AppendMe(c("chryso", "agi", "gale"))

fit = manova(log(as.matrix(chryso[,50:88])) ~ chryso$SEX + chryso$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercocebus_chrysogaster.txt", col.names = NA)

# Cercocebus galeritus
galeri = msrs[which(msrs$SPECIES == "galeritus"),]
#grupo irmao
agilis = msrs[which(msrs$SPECIES == "agilis"),]
agilis = agilis[which(agilis$GENUS == "Cercocebus"),]

galeri = rbind(galeri, agilis)

fit = manova(log(as.matrix(galeri[,50:88])) ~ galeri$SEX + galeri$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercocebus_galeritus.txt", col.names = NA)

# Cercocebus torquatus
torqua = msrs[which(msrs$SPECIES == "torquatus"),]

fit = lm(log(as.matrix(torqua[,50:88])) ~ torqua$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercocebus_torquatus.txt", col.names = NA)

# Cercopithecus albogularis
albogu = msrs[which(msrs$SPECIES == "albogularis"),]
albogu$SUBSPECIES[which(albogu$SUBSPECIES == "")] = "albogularis"

fit = manova(log(as.matrix(albogu[,50:88])) ~ albogu$SEX + albogu$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_albogularis.txt", col.names = NA)

# Cercopithecus ascanius
ascan = msrs[which(msrs$SPECIES == "ascanius"),]
ascan$SUBSPECIES[which(ascan$SUBSPECIES == "")] = "ascanius"

fit = manova(log(as.matrix(ascan[,50:88])) ~ ascan$SEX + ascan$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_ascanius.txt", col.names = NA)

# Cercopithecus campbelli
camp = msrs[which(msrs$SPECIES == "campbelli"),]

fit = lm(log(as.matrix(camp[,50:88])) ~ camp$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_campbelli.txt", col.names = NA)

# Cercopithecus cephus
cephus = msrs[which(msrs$SPECIES == "cephus"),]

fit = lm(log(as.matrix(cephus[,50:88])) ~ cephus$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_cephus.txt", col.names = NA)

# Cercopithecus denti
denti = msrs[which(msrs$SPECIES == "denti"),]

fit = lm(log(as.matrix(denti[,50:88])) ~ denti$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_denti.txt", col.names = NA)

# Cercopithecus diana
diana = msrs[which(msrs$SPECIES == "diana"),]

fit = lm(log(as.matrix(diana[,50:88])) ~ diana$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_diana.txt", col.names = NA)

# Cercopithecus doggetti
dog = msrs[which(msrs$SPECIES == "doggetti"),]
#grupo irmao
mitis = msrs[which(msrs$SPECIES == "mitis"),]
mitis$SUBSPECIES[which(mitis$SUBSPECIES == "")] = "mitis"
kan = msrs[which(msrs$SPECIES == "kandti"),]

dog = rbind(dog, mitis)
dog = rbind(dog, kan)

fit = manova(log(as.matrix(dog[,50:88])) ~ dog$SEX + dog$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_doggetti.txt", col.names = NA)

# Cercopithecus erythrogaster
ery = msrs[which(msrs$SPECIES == "erythrogaster"),]
#grupo irmao
peta = msrs[which(msrs$SPECIES == "petaurista"),]
peta$SUBSPECIES[which(peta$SUBSPECIES == "")] = "petaurista"

ery = rbind(ery, peta)

fit = manova(log(as.matrix(ery[,50:88])) ~ ery$SEX + ery$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_erythrogaster.txt", col.names = NA)

# Cercopithecus erythrotis
eryt = msrs[which(msrs$SPECIES == "erythrotis"),]
eryt$SUBSPECIES[which(eryt$SUBSPECIES == "")] = "erythrotis"
#eryt = eryt[-c(8,9),]
# grupo irmao
cephus = msrs[which(msrs$SPECIES == "cephus"),]

ascan = msrs[which(msrs$SPECIES == "ascanius"),]
ascan$SUBSPECIES[which(ascan$SUBSPECIES == "")] = "ascanius"

eryt = rbind(eryt, cephus)
eryt = rbind(eryt, ascan)

fit = manova(log(as.matrix(eryt[,50:88])) ~ eryt$SEX + eryt$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_erythrotis.txt", col.names = NA)

# Cercopithecus hamlyni
ham = msrs[which(msrs$SPECIES == "hamlyni"),]

fit = lm(log(as.matrix(ham[,50:88])) ~ ham$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_hamlyni.txt", col.names = NA)

# Cercopithecus kandti SINGLE MALE NO EFFECTS
kan = msrs[which(msrs$SPECIES == "kandti"),]
# grupo irmao
mitis = msrs[which(msrs$SPECIES == "mitis"),]
mitis$SUBSPECIES[which(mitis$SUBSPECIES == "")] = "mitis"

kan = rbind(kan, mitis)

fit = manova(log(as.matrix(kan[,50:88])) ~ kan$SEX + kan$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_kandti.txt", col.names = NA)

# Cercopithecus lhoesti
lhoes = msrs[which(msrs$SPECIES == "lhoesti"),]

fit = lm(log(as.matrix(lhoes[,50:88])) ~ lhoes$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_lhoesti.txt", col.names = NA)

# Cercopithecus lowei
lowei = msrs[which(msrs$SPECIES == "lowei"),]

fit = lm(log(as.matrix(lowei[,50:88])) ~ lowei$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_lowei.txt", col.names = NA)

# Cercopithecus mitis
mitis = msrs[which(msrs$SPECIES == "mitis"),]
mitis$SUBSPECIES[which(mitis$SUBSPECIES == "")] = "mitis"

fit = manova(log(as.matrix(mitis[,50:88])) ~ mitis$SEX + mitis$SUBSPECIES + 
               mitis$SEX:mitis$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_mitis.txt", col.names = NA)


# Cercopithecus mona
mona = msrs[which(msrs$SPECIES == "mona"),]

fit = lm(log(as.matrix(mona[,50:88])) ~ mona$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_mona.txt", col.names = NA)

# Cercopithecus neglectus
negl = msrs[which(msrs$SPECIES == "neglectus"),]

fit = lm(log(as.matrix(negl[,50:88])) ~ negl$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_neglectus.txt", col.names = NA)

# Cercopithecus nictitans
nict = msrs[which(msrs$SPECIES == "nictitans"),]
nict$SUBSPECIES[which(nict$SUBSPECIES == "")] = "nictitans"

fit = lm(log(as.matrix(nict[,50:88])) ~ nict$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_nictitans.txt", col.names = NA)

# Cercopithecus petaurista
peta = msrs[which(msrs$SPECIES == "petaurista"),]
peta$SUBSPECIES[which(peta$SUBSPECIES == "")] = "petaurista"

fit = lm(log(as.matrix(peta[,50:88])) ~ peta$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_petaurista.txt", col.names = NA)

# Cercopithecus pogonias
pogo = msrs[which(msrs$SPECIES == "pogonias"),]
pogo$SUBSPECIES[which(pogo$SUBSPECIES == "")] = "pogonias"

fit = lm(log(as.matrix(pogo[,50:88])) ~ pogo$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_pogonias.txt", col.names = NA)

# Cercopithecus preussi
preus = msrs[which(msrs$SPECIES == "preussi"),]

fit = lm(log(as.matrix(preus[,50:88])) ~ preus$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_preussi.txt", col.names = NA)

# Cercopithecus roloway
rolo = msrs[which(msrs$SPECIES == "roloway"),]
# grupo irmao
negl = msrs[which(msrs$SPECIES == "neglectus"),]
camp = msrs[which(msrs$SPECIES == "campbelli"),]
pogo = msrs[which(msrs$SPECIES == "pogonias"),]
pogo$SUBSPECIES[which(pogo$SUBSPECIES == "")] = "pogonias"
mona = msrs[which(msrs$SPECIES == "mona"),]
wolf = msrs[which(msrs$SPECIES == "wolfi"),]
denti = msrs[which(msrs$SPECIES == "denti"),]


rolo = AppendMe(c("rolo", "negl", "camp", "pogo", "mona", "wolf", "denti"))

fit = manova(log(as.matrix(rolo[,50:88])) ~ rolo$SEX + rolo$SPECIES +
               rolo$SEX:rolo$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_roloway.txt", col.names = NA)

# Cercopithecus wolfi
wolf = msrs[which(msrs$SPECIES == "wolfi"),]

fit = lm(log(as.matrix(wolf[,50:88])) ~ wolf$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Cercopithecus_wolfi.txt", col.names = NA)

# Chlorocebus aethiops
aeth = msrs[which(msrs$SPECIES == "aethiops"),]
aeth$SUBSPECIES[which(aeth$SUBSPECIES == "0")] = "aethiops"

fit = manova(log(as.matrix(aeth[,50:88])) ~ aeth$SEX + aeth$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Chlorocebus_aethiops.txt", col.names = NA)

# Chlorocebus cynosuros
cyno = msrs[which(msrs$SPECIES == "cynosuros"),]
# grupo irmao
pyge = msrs[which(msrs$SPECIES == "pygerythrus"),]
pyge$SUBSPECIES[which(pyge$SUBSPECIES == "")] = "pygerythrus"

cyno = AppendMe(c("cyno", "pyge"))

fit = manova(as.matrix(cyno[,50:88]) ~ cyno$SEX + cyno$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Chlorocebus_cynosuros.txt", col.names = NA)

# Chlorocebus djamdjamensis
djam = msrs[which(msrs$SPECIES == "djamdjamensis"),]
# grupo irmao
cyno = msrs[which(msrs$SPECIES == "cynosuros"),]
pyge = msrs[which(msrs$SPECIES == "pygerythrus"),]
pyge$SUBSPECIES[which(pyge$SUBSPECIES == "")] = "pygerythrus"
aeth = msrs[which(msrs$SPECIES == "aethiops"),]
aeth$SUBSPECIES[which(aeth$SUBSPECIES == "0")] = "aethiops"
tan = msrs[which(msrs$SPECIES == "tantalus"),]
tan$SUBSPECIES[which(tan$SUBSPECIES == "")] = "tantalus"

djam = AppendMe(c("djam", "cyno", "pyge", "aeth", "tan"))

fit = manova(log(as.matrix(djam[,50:88])) ~ djam$SEX + djam$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Chlorocebus_djamdjamensis.txt", col.names = NA)

# Chlorocebus pygerythrus
pyge = msrs[which(msrs$SPECIES == "pygerythrus"),]
pyge$SUBSPECIES[which(pyge$SUBSPECIES == "")] = "pygerythrus"

fit = manova(log(as.matrix(pyge[,50:88])) ~ pyge$SEX + pyge$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Chlorocebus_pygerythrus.txt", col.names = NA)

# Chlorocebus sabaeus
saba = msrs[which(msrs$SPECIES == "sabaeus"),]
  
fit = lm(log(as.matrix(saba[,50:88])) ~ saba$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Chlorocebus_sabaeus.txt", col.names = NA)

# Chlorocebus tantalus
tan = msrs[which(msrs$SPECIES == "tantalus"),]
tan$SUBSPECIES[which(tan$SUBSPECIES == "")] = "tantalus"

fit = lm(log(as.matrix(tan[,50:88])) ~ tan$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Chlorocebus_tantalus.txt", col.names = NA)

# Colobus angolensis
ango = msrs[which(msrs$SPECIES == "angolensis"),]
ango$SUBSPECIES[which(ango$SUBSPECIES == "")] = "angolensis"

fit = manova(log(as.matrix(ango[,50:88])) ~ ango$SEX + ango$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Colobus_angolensis.txt", col.names = NA)

# Colobus guereza
guere = msrs[which(msrs$SPECIES == "guereza"),]
guere$SUBSPECIES[which(guere$SUBSPECIES == "")] = "guereza"

fit = manova(log(as.matrix(guere[,50:88])) ~ guere$SEX + guere$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Colobus_guereza.txt", col.names = NA)

# Colobus polykomos
poly = msrs[which(msrs$SPECIES == "polykomos"),]
# grupo irmao
vele = msrs[which(msrs$SPECIES == "vellerosus"),]

poly = rbind(poly, vele)

fit = manova(log(as.matrix(poly[,50:88])) ~ poly$SEX + poly$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Colobus_polykomos.txt", col.names = NA)

# Colobus satanas
sata = msrs[which(msrs$SPECIES == "satanas"),]
sata$SUBSPECIES[which(sata$SUBSPECIES == "")] = "satanas"
sata$SUBSPECIES[which(sata$SUBSPECIES == "0")] = "satanas"
# grupo irmao
ango = msrs[which(msrs$SPECIES == "angolensis"),]
ango$SUBSPECIES[which(ango$SUBSPECIES == "")] = "angolensis"

sata = rbind(sata, ango)

fit = manova(log(as.matrix(sata[,50:88])) ~ sata$SEX + sata$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Colobus_satanas.txt", col.names = NA)


# Colobus vellerosus
vele = msrs[which(msrs$SPECIES == "vellerosus"),]

fit = lm(log(as.matrix(vele[,50:88])) ~ vele$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Colobus_vellerosus.txt", col.names = NA)

# Erythrocebus patas
patas = msrs[which(msrs$SPECIES == "patas"),]

fit = lm(log(as.matrix(patas[,50:88])) ~ patas$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Erythrocebus_patas.txt", col.names = NA)

# Gorilla beringei
gor = msrs[which(msrs$SPECIES == "beringei"),]
gor$SUBSPECIES[which(gor$SUBSPECIES == "")] = "beringei"

fit = manova(log(as.matrix(gor[,50:88])) ~ gor$SEX + gor$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Gorilla_beringei.txt", col.names = NA)

# Gorilla gorilla
gori = msrs[which(msrs$SPECIES == "gorilla"),]

fit = lm(log(as.matrix(gori[,50:88])) ~ gori$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Gorilla_gorilla.txt", col.names = NA)

# Homo sapiens
homo = msrs[which(msrs$SPECIES == "sapiens"),]

fit = lm(log(as.matrix(homo[,50:88])) ~ homo$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Homo_sapiens.txt", col.names = NA)

# Hylobates agilis
hylo = msrs[which(msrs$GENUS == "Hylobates"),]
hylo = hylo[which(hylo$SPECIES == "agilis"),]

hylo = hylo[-which(hylo$SEX == "sexo"),]
hylo = hylo[-which(hylo$SEX == "?male"),]
            
fit = lm(log(as.matrix(hylo[,50:88])) ~ hylo$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Hylobates_agilis.txt", col.names = NA)

# Hylobates albibarbis
albi = msrs[which(msrs$SPECIES == "albibarbis"),]
# grupo irmao
hylo = msrs[which(msrs$GENUS == "Hylobates"),]
hylo = hylo[which(hylo$SPECIES == "agilis"),]

hylo = hylo[-which(hylo$SEX == "sexo"),]
hylo = hylo[-which(hylo$SEX == "?male"),]

albi = rbind(albi, hylo)

fit = manova(log(as.matrix(albi[,50:88])) ~ albi$SEX + albi$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Hylobates_albibarbis.txt", col.names = NA)

# Hylobates klossi
klos = msrs[which(msrs$SPECIES == "klossii"),]
# grupo irmao
mol = msrs[which(msrs$SPECIES == "moloch"),]
mol = mol[-which(mol$SEX == "sexo"),]
mol = mol[-which(mol$SEX == "?female"),]
albi = msrs[which(msrs$SPECIES == "albibarbis"),]
hylo = msrs[which(msrs$GENUS == "Hylobates"),]
hylo = hylo[which(hylo$SPECIES == "agilis"),]

hylo = hylo[-which(hylo$SEX == "sexo"),]
hylo = hylo[-which(hylo$SEX == "?male"),]

klos = AppendMe(c("klos", "mol", "albi", "hylo"))

fit = manova(log(as.matrix(klos[,50:88])) ~ klos$SEX + klos$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Hylobates_klossii.txt", col.names = NA)

# Hylobates lar
lar = msrs[which(msrs$SPECIES == "lar"),]

fit = manova(log(as.matrix(lar[,50:88])) ~ lar$SEX + lar$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Hylobates_lar.txt", col.names = NA)

# Hylobates moloch
mol = msrs[which(msrs$SPECIES == "moloch"),]

mol = mol[-which(mol$SEX == "sexo"),]
mol = mol[-which(mol$SEX == "?female"),]
# grupo irmao
klos = msrs[which(msrs$SPECIES == "klossii"),]
albi = msrs[which(msrs$SPECIES == "albibarbis"),]
hylo = msrs[which(msrs$GENUS == "Hylobates"),]
hylo = hylo[which(hylo$SPECIES == "agilis"),]

hylo = hylo[-which(hylo$SEX == "sexo"),]
hylo = hylo[-which(hylo$SEX == "?male"),]

mol = AppendMe(c("mol", "klos", "albi", "hylo"))

fit = manova(log(as.matrix(mol[,50:88])) ~ mol$SEX + mol$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Hylobates_moloch.txt", col.names = NA)

# Hylobates muelleri
mue = msrs[which(msrs$SPECIES == "muelleri"),]

mue$SUBSPECIES[which(mue$SUBSPECIES == "")] = "muelleri"
mue$SUBSPECIES[which(mue$SUBSPECIES == "0")] = "muelleri"

mue = mue[-which(mue$SEX == "sexo"),]
mue = mue[-which(mue$SEX == "?male"),]
mue = mue[-which(mue$SEX == "0"),]

fit = lm(log(as.matrix(mue[,50:88])) ~ mue$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Hylobates_muelleri.txt", col.names = NA)

# Hylobateds pileatus
pile = msrs[which(msrs$GENUS == "Hylobates"),]
#pile = pile[which(pile$SPECIES == "pileatus"),]
# grupo irmao usar todo o genero Hylobates

fit = manova(log(as.matrix(pile[,50:88])) ~ pile$SEX + pile$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Hylobates_pileatus.txt", col.names = NA)

# Kasi johnii
john = msrs[which(msrs$SPECIES == "johnii"),]
# grupo irmao
ente = msrs[which(msrs$SPECIES == "entellus"),]
hec = msrs[which(msrs$SPECIES == "hector"),]

john = AppendMe(c("john", "ente", "hec"))

fit = manova(log(as.matrix(john[,50:88])) ~ john$SEX + john$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Kasi_johnii.txt", col.names = NA)

# Kasi vetulus
vetu = msrs[which(msrs$SPECIES == "vetulus"),]
vetu$SUBSPECIES[which(vetu$SUBSPECIES == "")] = "vetulus"
# grupo irmao
pri = msrs[which(msrs$SPECIES == "priam"),]

vetu = rbind(vetu, pri)

fit = manova(log(as.matrix(vetu[,50:88])) ~ vetu$SEX + vetu$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Kasi_vetulus.txt", col.names = NA)

# Lophocebus albigena
lopho = msrs[which(msrs$GENUS == "Lophocebus"),]
lopho = lopho[which(lopho$SPECIES == "albigena"),]

lopho$SUBSPECIES[which(lopho$SUBSPECIES == "")] = "albigena"

fit = manova(log(as.matrix(lopho[,50:88])) ~ lopho$SEX + lopho$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Lophocebus_albigena.txt", col.names = NA)

# Lophocebus aterrimus
lopho = msrs[which(msrs$GENUS == "Lophocebus"),]
lopho = lopho[which(lopho$SPECIES == "aterrimus"),]

fit = lm(log(as.matrix(lopho[,50:88])) ~ lopho$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Lophocebus_aterrimus.txt", col.names = NA)

# Lophocebus opdenboschi
opden = msrs[which(msrs$GENUS == "Lophocebus"),]
opden = opden[which(opden$SPECIES == "opdenboschi"),]
# grupo irmao
lopho = msrs[which(msrs$GENUS == "Lophocebus"),]
lopho = lopho[which(lopho$SPECIES == "aterrimus"),]

opden = rbind(opden, lopho)

fit = manova(log(as.matrix(opden[,50:88])) ~ opden$SEX + opden$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Lophocebus_opdenboschi.txt", col.names = NA)

# Macaca arctoides
arcto = msrs[which(msrs$SPECIES == "arctoides"),]

fit = lm(log(as.matrix(arcto[,50:88])) ~ arcto$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_arctoides.txt", col.names = NA)

# Macaca assamensis
assa = msrs[which(msrs$SPECIES == "assamensis"),]
assa$SUBSPECIES[which(assa$SUBSPECIES == "")] = "assamensis"
assa$SUBSPECIES[which(assa$SUBSPECIES == "0")] = "assamensis"

fit = lm(log(as.matrix(assa[,50:88])) ~ assa$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_assamensis.txt", col.names = NA)

# Macaca cyclops
cyclo = msrs[which(msrs$SPECIES == "cyclopis"),]

fit = lm(log(as.matrix(cyclo[,50:88])) ~ cyclo$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_cyclopis.txt", col.names = NA)

# Macaca fascicularis
fasci = msrs[which(msrs$SPECIES == "fascicularis"),]
fasci$SUBSPECIES[which(fasci$SUBSPECIES == "")] = "fascicularis"

fit = manova(log(as.matrix(fasci[,50:88])) ~ fasci$SEX + fasci$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_fascicularis.txt", col.names = NA)

# Macaca fuscata
fusca = msrs[which(msrs$SPECIES == "fuscata"),]

fit = lm(log(as.matrix(fusca[,50:88])) ~ fusca$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_fuscata.txt", col.names = NA)

# Macaca hecki
heck = msrs[which(msrs$SPECIES == "hecki"),]
# grupo irmao
maura = msrs[which(msrs$SPECIES == "maura"),]
och = msrs[which(msrs$SPECIES == "ochreata"),]

heck = AppendMe(c("heck", "maura", "och"))

fit = manova(log(as.matrix(heck[,50:88])) ~ heck$SEX + heck$SPECIES +
               heck$SEX:heck$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_hecki.txt", col.names = NA)

# Macaca leonina
leo = msrs[which(msrs$SPECIES == "leonina"),]
# grupo irmao
page = msrs[which(msrs$SPECIES == "pagensis"),]
sile = msrs[which(msrs$SPECIES == "silenus"),]

leo = AppendMe(c("leo", "page", "sile"))

fit = manova(log(as.matrix(leo[,50:88])) ~ leo$SEX + leo$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_leonina.txt", col.names = NA)

# Macaca maura
maura = msrs[which(msrs$SPECIES == "maura"),]
# grupo irmao
tonk = msrs[which(msrs$SPECIES == "tonkeana"),]
nigra = msrs[which(msrs$SPECIES == "nigra"),]

maura = AppendMe(c("maura", "tonk", "nigra"))

fit = manova(as.matrix(log(maura[,50:88])) ~ maura$SEX + maura$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_maura.txt", col.names = NA)

# Macaca mulatta
mula = msrs[which(msrs$SPECIES == "mulatta"),]

fit = lm(log(as.matrix(mula[,50:88])) ~ mula$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_mulatta.txt", col.names = NA)

# Macaca nemestrina
neme = msrs[which(msrs$SPECIES == "nemestrina"),]

fit = lm(log(as.matrix(neme[,50:88])) ~ neme$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_nemestrina.txt", col.names = NA)

# Macaca nigra
nigra = msrs[which(msrs$SPECIES == "nigra"),]

fit = lm(log(as.matrix(nigra[,50:88])) ~ nigra$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_nigra.txt", col.names = NA)

# Macaca ochreata
och = msrs[which(msrs$SPECIES == "ochreata"),]
# grupo irmao
maura = msrs[which(msrs$SPECIES == "maura"),]

och = rbind(och, maura)

fit = manova(log(as.matrix(och[,50:88])) ~ och$SEX + och$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_ochreata.txt", col.names = NA)

# Macaca pagensis
page = msrs[which(msrs$SPECIES == "pagensis"),]
# grupo irmao
leo = msrs[which(msrs$SPECIES == "leonina"),]
sile = msrs[which(msrs$SPECIES == "silenus"),]

page = AppendMe(c("page", "leo", "sile"))

fit = manova(log(as.matrix(page[,50:88])) ~ page$SEX + page$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_pagensis.txt", col.names = NA)

# Macaca radiata
radi = msrs[which(msrs$SPECIES == "radiata"),]

fit = lm(log(as.matrix(radi[,50:88])) ~ radi$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_radiata.txt", col.names = NA)

# Macaca silenus
sile = msrs[which(msrs$SPECIES == "silenus"),]
# grupo irmao
page = msrs[which(msrs$SPECIES == "pagensis"),]
leo = msrs[which(msrs$SPECIES == "leonina"),]

sile = AppendMe(c("sile", "page", "leo"))

fit = manova(log(as.matrix(sile[,50:88])) ~ sile$SEX + sile$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_silenus.txt", col.names = NA)

# Macaca sinica
sini = msrs[which(msrs$SPECIES == "sinica"),]

fit = lm(log(as.matrix(sini[,50:88])) ~ sini$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_sinica.txt", col.names = NA)

# Macaca sylvanus
sylva = msrs[which(msrs$SPECIES == "sylvanus"),]

fit = lm(log(as.matrix(sylva[,50:88])) ~ sylva$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_sylvanus.txt", col.names = NA)

# Macaca thibetana
tybe = msrs[which(msrs$SPECIES == "thibetana"),]
# grupo irmao
assa = msrs[which(msrs$SPECIES == "assamensis"),]
assa$SUBSPECIES[which(assa$SUBSPECIES == "")] = "assamensis"
assa$SUBSPECIES[which(assa$SUBSPECIES == "0")] = "assamensis"

tybe = rbind(tybe, assa)

fit = manova(log(as.matrix(tybe[,50:88])) ~ tybe$SEX + tybe$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_thibetana.txt", col.names = NA)

# Macaca tonkeana
tonk = msrs[which(msrs$SPECIES == "tonkeana"),]
# grupo irmao
nigra = msrs[which(msrs$SPECIES == "nigra"),]

tonk = rbind(tonk, nigra)

fit = manova(log(as.matrix(tonk[,50:88])) ~ tonk$SEX + tonk$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Macaca_tonkeana.txt", col.names = NA)

# Mandrillus leucophaeus
mand = msrs[which(msrs$SPECIES == "leucophaeus"),]

fit = lm(log(as.matrix(mand[,50:88])) ~ mand$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Mandrillus_leucophaeus.txt", col.names = NA)

# Mandrillus sphinx
sph = msrs[which(msrs$SPECIES == "sphinx"),]

fit = lm(log(as.matrix(sph[,50:88])) ~ sph$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Mandrillus_sphinx.txt", col.names = NA)

# Miopithecus ogouensis
mio = msrs[which(msrs$SPECIES == "ogouensis"),]
# grupo irmao
tala = msrs[which(msrs$SPECIES == "talapoin"),]

mio = AppendMe(c("mio", "tala"))

fit = lm(log(as.matrix(mio[,50:88])) ~ mio$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Miopithecus_ogouensis.txt", col.names = NA)

# Miopithecus talapoin
mio = msrs[which(msrs$SPECIES == "talapoin"),]

fit = lm(log(as.matrix(mio[,50:88])) ~ mio$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Miopithecus_talapoin.txt", col.names = NA)

# Nasalis larvatus
nasa = msrs[which(msrs$SPECIES == "larvatus"),]

fit = lm(log(as.matrix(nasa[,50:88])) ~ nasa$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Nasalis_larvatus.txt", col.names = NA)

# Nomascus concolor
conco = msrs[which(msrs$SPECIES == "concolor"),]

fit = lm(log(as.matrix(conco[,50:88])) ~ conco$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Nomascus_concolor.txt", col.names = NA)

# Nomascus gabriellae
gab = msrs[which(msrs$SPECIES == "gabriellae"),]
# grupo irmao
leuco = msrs[which(msrs$SPECIES == "leucogenys"),]
leuco = leuco[-which(leuco$SEX == "?female"),]

gab = rbind(gab, leuco)

fit = manova(log(as.matrix(gab[,50:88])) ~ gab$SEX + gab$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Nomascus_gabriellae.txt", col.names = NA)

# Nomascus leucogenys
leuco = msrs[which(msrs$SPECIES == "leucogenys"),]
leuco = leuco[-which(leuco$SEX == "?female"),]
# grupo irmao
gab = msrs[which(msrs$SPECIES == "gabriellae"),]

leuco = rbind(leuco, gab)

fit = manova(log(as.matrix(leuco[,50:88])) ~ leuco$SEX + leuco$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Nomascus_leucogenys.txt", col.names = NA)

# Pan paniscus
pan = msrs[which(msrs$SPECIES == "paniscus"),]

fit = lm(log(as.matrix(pan[,50:88])) ~ pan$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Pan_paniscus.txt", col.names = NA)

# Pan troglodytes
tro = msrs[which(msrs$SPECIES == "troglodytes"),]
tro$SUBSPECIES[which(tro$SUBSPECIES == "")] = "troglodytes"

fit = manova(log(as.matrix(tro[,50:88])) ~ tro$SEX + tro$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Pan_troglodytes.txt", col.names = NA)

# Papio anubis
anu = msrs[which(msrs$SPECIES == "anubis"),]

fit = lm(log(as.matrix(anu[,50:88])) ~ anu$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Papio_anubis.txt", col.names = NA)

# Papio cynocephalus
cyno = msrs[which(msrs$SPECIES == "cynocephalus"),]
cyno$SUBSPECIES[which(cyno$SUBSPECIES == "")] = "cynocephalus"

fit = manova(log(as.matrix(cyno[,50:88])) ~ cyno$SEX + cyno$SUBSPECIES +
           cyno$SEX:cyno$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Papio_cynocephalus.txt", col.names = NA)

# Papio hamadryas
hama = msrs[which(msrs$SPECIES == "hamadryas"),]
hama$SUBSPECIES[which(hama$SUBSPECIES == "")] = "hamadryas"
hama$SUBSPECIES[which(hama$SUBSPECIES == "0")] = "hamadryas"

fit = lm(log(as.matrix(hama[,50:88])) ~ hama$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Papio_hamadryas.txt", col.names = NA)

# Papio papio
papio = msrs[which(msrs$SPECIES == "papio"),]
# grupo irmao
anu = msrs[which(msrs$SPECIES == "anubis"),]
hama = msrs[which(msrs$SPECIES == "hamadryas"),]
hama$SUBSPECIES[which(hama$SUBSPECIES == "")] = "hamadryas"
hama$SUBSPECIES[which(hama$SUBSPECIES == "0")] = "hamadryas"

papio = AppendMe(c("papio", "anu", "hama"))

fit = manova(log(as.matrix(papio[,50:88])) ~ papio$SEX + papio$SPECIES +
               papio$SEX:papio$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Papio_papio.txt", col.names = NA)

# Papio ursinus
ursi = msrs[which(msrs$SPECIES == "ursinus"),]
ursi$SUBSPECIES[which(ursi$SUBSPECIES == "")] = "ursinus"

fit = lm(log(as.matrix(ursi[,50:88])) ~ ursi$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Papio_ursinus.txt", col.names = NA)

# Piliocolobus badius
badi = msrs[which(msrs$SPECIES == "badius"),]
badi$SUBSPECIES[which(badi$SUBSPECIES == "")] = "badius"

fit = manova(log(as.matrix(badi[,50:88])) ~ badi$SEX + badi$SUBSPECIES + badi$SEX:badi$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_badius.txt", col.names = NA)

# Piliocolobus foai
foai = msrs[which(msrs$SPECIES == "foai"),]
foai$SUBSPECIES[which(foai$SUBSPECIES == "")] = "foai"

fit = manova(log(as.matrix(foai[,50:88])) ~ foai$SEX + foai$SUBSPECIES + foai$SEX:foai$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_foai.txt", col.names = NA)

# Piliocolobus kirkii
kirki = msrs[which(msrs$SPECIES == "kirkii"),]
# grupo irmao
foai = msrs[which(msrs$SPECIES == "foai"),]
foai$SUBSPECIES[which(foai$SUBSPECIES == "")] = "foai"
tep = msrs[which(msrs$SPECIES == "tephrosceles"),]
rufo = msrs[which(msrs$SPECIES == "rufomitratus"),]

kirki = AppendMe(c("kirki", "foai", "tep", "rufo"))

fit = manova(log(as.matrix(kirki[,50:88])) ~ kirki$SEX + kirki$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_kirkii.txt", col.names = NA)

# Piliocolobus pennantii
pen = msrs[which(msrs$SPECIES == "pennantii"),]
pen$SUBSPECIES[which(pen$SUBSPECIES == "")] = "pennantii"

fit = manova(log(as.matrix(pen[,50:88])) ~ pen$SEX + pen$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_pennantii.txt", col.names = NA)

# Piliocolobus preussi
preu = msrs[which(msrs$SPECIES == "preussi"),]
preu = preu[which(preu$GENUS == "Piliocolobus"),]
# grupo irmao
tho = msrs[which(msrs$SPECIES == "tholloni"),]

preu = rbind(preu, tho)

fit = manova(log(as.matrix(preu[,50:88])) ~ preu$SEX + preu$SPECIES +
               preu$SEX:preu$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_preussi.txt", col.names = NA)

# Piliocolobus rufomitratus
rufo = msrs[which(msrs$SPECIES == "rufomitratus"),]
# grupo irmao
tep = msrs[which(msrs$SPECIES == "tephrosceles"),]

rufo = rbind(rufo, tep)

fit = manova(log(as.matrix(rufo[,50:88])) ~ rufo$SEX + rufo$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_rufomitratus.txt", col.names = NA)

# Piliocolobus tephrosceles
tep = msrs[which(msrs$SPECIES == "tephrosceles"),]
# grupo irmao
rufo = msrs[which(msrs$SPECIES == "rufomitratus"),]

tep = rbind(tep, rufo)

fit = manova(log(as.matrix(tep[,50:88])) ~ tep$SEX + tep$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_tephrosceles.txt", col.names = NA)

# Piliocolobus tholloni
tho = msrs[which(msrs$SPECIES == "tholloni"),]

fit = lm(log(as.matrix(tho[,50:88])) ~ tho$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Piliocolobus_tholloni.txt", col.names = NA)

# Pongo abelii
abe = msrs[which(msrs$SPECIES == "abelii"),]
# grupo irmao
pyg = msrs[which(msrs$SPECIES == "pygmaeus"),]

abe = rbind(abe, pyg)

fit = manova(log(as.matrix(abe[,50:88])) ~ abe$SEX + abe$SPECIES +
               abe$SEX:abe$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Pongo_abelii.txt", col.names = NA)

# Pongo pygmaeus
pyg = msrs[which(msrs$SPECIES == "pygmaeus"),]

fit = lm(log(as.matrix(pyg[,50:88])) ~ pyg$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Pongo_pygmaeus.txt", col.names = NA)

# Presbytis chrysomelas
chryso = msrs[which(msrs$SPECIES == "chrysomelas"),]
chryso$SUBSPECIES[which(chryso$SUBSPECIES == "")] = "chrysomelas"

fit = lm(log(as.matrix(chryso[,50:88])) ~ chryso$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_chrysomelas.txt", col.names = NA)

# Presbytis comata
coma = msrs[which(msrs$SPECIES == "comata"),]

fit = lm(log(as.matrix(coma[,50:88])) ~ coma$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_comata.txt", col.names = NA)

# Presbytis femoralis
femo = msrs[which(msrs$SPECIES == "femoralis"),]
femo$SUBSPECIES[which(femo$SUBSPECIES == "")] = "femoralis"

fit = lm(log(as.matrix(femo[,50:88])) ~ femo$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_femoralis.txt", col.names = NA)

# Presbytis frontata
fron = msrs[which(msrs$SPECIES == "frontata"),]
# grupo irmao
chryso = msrs[which(msrs$SPECIES == "chrysomelas"),]
chryso$SUBSPECIES[which(chryso$SUBSPECIES == "")] = "chrysomelas"

fron = rbind(fron, chryso)

fit = lm(log(as.matrix(fron[,50:88])) ~ fron$SEX + fron$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_frontata.txt", col.names = NA)

# Presbytis hosei
hosei = msrs[which(msrs$SPECIES == "hosei"),]
# grupo irmao
fron = msrs[which(msrs$SPECIES == "frontata"),]
chryso = msrs[which(msrs$SPECIES == "chrysomelas"),]
chryso$SUBSPECIES[which(chryso$SUBSPECIES == "")] = "chrysomelas"

hosei = AppendMe(c("hosei", "fron", "chryso"))

fit = manova(log(as.matrix(hosei[,50:88])) ~ hosei$SEX + hosei$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_hosei.txt", col.names = NA)

# Presbytis melalophos
mela = msrs[which(msrs$SPECIES == "melalophos"),]

fit = lm(log(as.matrix(mela[,50:88])) ~ mela$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_melalophos.txt", col.names = NA)

# Presbytis natunae
natu = msrs[which(msrs$SPECIES == "natunae"),]

fit = lm(log(as.matrix(natu[,50:88])) ~ natu$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_natunae.txt", col.names = NA)

# Presbytis potenziani
pote = msrs[which(msrs$SPECIES == "potenziani"),]
# grupo irmao
coma = msrs[which(msrs$SPECIES == "comata"),]
rubi = msrs[which(msrs$SPECIES == "rubicunda"),]
rubi$SUBSPECIES[which(rubi$SUBSPECIES == "")] = "rubicunda"
mela = msrs[which(msrs$SPECIES == "melalophos"),]

pote = AppendMe(c("pote", "coma", "rubi", "mela"))

fit = manova(log(as.matrix(pote[,50:88])) ~ pote$SEX + pote$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_potenziani.txt", col.names = NA)

# Presbytis rubicunda
rubi = msrs[which(msrs$SPECIES == "rubicunda"),]
rubi$SUBSPECIES[which(rubi$SUBSPECIES == "")] = "rubicunda"

fit = manova(log(as.matrix(rubi[,50:88])) ~ rubi$SEX + rubi$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_rubicunda.txt", col.names = NA)

# Presbytis siamensis
sia = msrs[which(msrs$SPECIES == "siamensis"),]
sia$SUBSPECIES[which(sia$SUBSPECIES == "")] = "siamensis"

fit = lm(log(as.matrix(sia[,50:88])) ~ sia$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_siamensis.txt", col.names = NA)

# Presbytis thomasi
tho = msrs[which(msrs$SPECIES == "thomasi"),]
# grupo irmao
hosei = msrs[which(msrs$SPECIES == "hosei"),]
fron = msrs[which(msrs$SPECIES == "frontata"),]
chryso = msrs[which(msrs$SPECIES == "chrysomelas"),]
chryso$SUBSPECIES[which(chryso$SUBSPECIES == "")] = "chrysomelas"

tho = AppendMe(c("tho", "hosei", "fron", "chryso"))

fit = manova(log(as.matrix(tho[,50:88])) ~ tho$SEX + tho$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Presbytis_thomasi.txt", col.names = NA)

# Procolobus verus
verus = msrs[which(msrs$SPECIES == "verus"),]

fit = lm(log(as.matrix(verus[,50:88])) ~ verus$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Procolobus_verus.txt", col.names = NA)

# Pygathrix nemaeus
nema = msrs[which(msrs$SPECIES == "nemaeus"),]
# grupo irmao
nigri = msrs[which(msrs$SPECIES == "nigripes"),]

nema = rbind(nema, nigri)

fit = manova(log(as.matrix(nema[,50:88])) ~ nema$SEX + nema$SPECIES +
               nema$SEX:nema$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Pygathrix_nemaeus.txt", col.names = NA)

# Pygathrix nigripes
nigri = msrs[which(msrs$SPECIES == "nigripes"),]
# grupo irmao
nema = msrs[which(msrs$SPECIES == "nemaeus"),]

nigri = rbind(nigri, nema)

fit = manova(log(as.matrix(nigri[,50:88])) ~ nigri$SEX + nigri$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Pygathrix_nigripes.txt", col.names = NA)

# Rhinopithecus avunculus
avu = msrs[which(msrs$SPECIES == "avunculus"),]
# grupo irmao
roxe = msrs[which(msrs$SPECIES == "roxellana"),]
con = msrs[which(msrs$SPECIES == "concolor"),]
nasa = msrs[which(msrs$SPECIES == "larvatus"),]

avu = AppendMe(c("avu", "roxe", "con", "nasa"))

fit = manova(log(as.matrix(avu[,50:88])) ~ avu$SEX + avu$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Rhinopithecus_avunculus.txt", col.names = NA)

# Rhinopithecus roxellana
roxe = msrs[which(msrs$SPECIES == "roxellana"),]
# grupo irmao
avu = msrs[which(msrs$SPECIES == "avunculus"),]
con = msrs[which(msrs$SPECIES == "concolor"),]
nasa = msrs[which(msrs$SPECIES == "larvatus"),]

roxe = AppendMe(c("avu", "con", "nasa"))

fit = manova(log(as.matrix(roxe[,50:88])) ~ roxe$SEX + roxe$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Rhinopithecus_roxellana.txt", col.names = NA)

# Semnopithecus ajax
ajax = msrs[which(msrs$SPECIES == "ajax"),]

fit = lm(log(as.matrix(ajax[,50:88])) ~ ajax$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Semnopithecus_ajax.txt", col.names = NA)

# Semnopithecus dussumieri
dussu = msrs[which(msrs$SPECIES == "dussumieri"),]

fit = lm(log(as.matrix(dussu[,50:88])) ~ dussu$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Semnopithecus_dussumieri.txt", col.names = NA)

# Semnopithecus entellus
ente = msrs[which(msrs$SPECIES == "entellus"),]

fit = lm(log(as.matrix(ente[,50:88])) ~ ente$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Semnopithecus_entellus.txt", col.names = NA)

# Semnopithecus hector
hec = msrs[which(msrs$SPECIES == "hector"),]
# grupo irmao
ente = msrs[which(msrs$SPECIES == "entellus"),]

hec = rbind(hec, ente)

fit = manova(log(as.matrix(hec[,50:88])) ~ hec$SEX + hec$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Semnopithecus_hector.txt", col.names = NA)

# Semnopithecus priam
pri = msrs[which(msrs$SPECIES == "priam"),]

fit = lm(log(as.matrix(pri[,50:88])) ~ pri$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Semnopithecus_priam.txt", col.names = NA)

# Simias concolor
con = msrs[which(msrs$SPECIES == "concolor"),]

fit = lm(log(as.matrix(con[,50:88])) ~ con$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Simias_concolor.txt", col.names = NA)

# Symphalangus syndactylus
synda = msrs[which(msrs$SPECIES == "syndactylus"),]
synda = synda[-which(synda$SEX == "sexo"),]
synda = synda[-which(synda$SEX == "0"),]

fit = lm(log(as.matrix(synda[,50:88])) ~ synda$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Symphalangus_syndactylus.txt", col.names = NA)

# Theropithecus gelada
gela = msrs[which(msrs$SPECIES == "gelada"),]
gela$SUBSPECIES[which(gela$SUBSPECIES == "")] = "gelada"

fit = manova(log(as.matrix(gela[,50:88])) ~ gela$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Theropithecus_gelada.txt", col.names = NA)

# Trachypithecus auratus
aura = msrs[which(msrs$SPECIES == "auratus"),]

fit = lm(log(as.matrix(aura[,50:88])) ~ aura$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_auratus.txt", col.names = NA)

# Trachypithecus cristatus
crist = msrs[which(msrs$SPECIES == "cristatus"),]

fit = lm(log(as.matrix(crist[,50:88])) ~ crist$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_cristatus.txt", col.names = NA)

# Trachypithecus francoisi
fran = msrs[which(msrs$SPECIES == "francoisi"),]

fit = lm(log(as.matrix(fran[,50:88])) ~ fran$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_francoisi.txt", col.names = NA)

# Trachypithecus germaini
ger = msrs[which(msrs$SPECIES == "germaini"),]
# grupo irmao
aura = msrs[which(msrs$SPECIES == "auratus"),]

ger = rbind(ger, aura)

fit = manova(log(as.matrix(ger[,50:88])) ~ ger$SEX + ger$SPECIES +
               ger$SEX:ger$SPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_germaini.txt", col.names = NA)

# Trachypithecus obscurus
obsc = msrs[which(msrs$SPECIES == "obscurus"),]

fit = manova(log(as.matrix(obsc[,50:88])) ~ obsc$SEX + obsc$SUBSPECIES)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_obscurus.txt", col.names = NA)

# Trachypithecus phayrei
pha = msrs[which(msrs$SPECIES == "phayrei"),]
pha$SUBSPECIES[which(pha$SUBSPECIES == "")] = "phayrei"

fit = lm(log(as.matrix(pha[,50:88])) ~ pha$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_phayrei.txt", col.names = NA)

# Trachypithecus pileatus
pile = msrs[which(msrs$SPECIES == "pileatus"),]
pile = pile[which(pile$GENUS == "Trachypithecus"),]
pile$SUBSPECIES[which(pile$SUBSPECIES == "")] = "pileatus"

fit = lm(log(as.matrix(pile[,50:88])) ~ pile$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_pileatus.txt", col.names = NA)

# Trachypithecus shortridgei
sho = msrs[which(msrs$SPECIES == "shortridgei"),]

fit = lm(log(as.matrix(sho[,50:88])) ~ sho$SEX)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/log_vcv/Trachypithecus_shortridgei.txt", col.names = NA)
