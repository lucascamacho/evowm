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
msrs = read.csv(file = "~/Dropbox/Doc/Data/primates_measures/medidas_platyrrhini.csv", dec = ",", sep = ",")

# new column with species names
species = vector()
for(i in 1:nrow(msrs)){
  if(is.na(msrs$SUB.[i])){
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], sep="_")}
  else{
    species[i] = paste(msrs$GENUS.[i], msrs$SPECIES.[i], msrs$SUB.[i], sep="_")
  }
}

msrs = cbind(species, msrs)
species = unique(species)

index = which(apply(msrs[,24:62], 1, anyNA) == TRUE)
msrs = msrs[-index,]

###############################################################################################################

# Alouatta belzebull
belze = msrs[which(msrs$SPECIES. == "belzebul"),]

fit = manova(as.matrix(belze[,24:62]) ~ belze$SEX4. + belze$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Alouatta_belzebul.txt", col.names = NA)

# Alouatta caraya
caray = msrs[which(msrs$SPECIES. == "caraya"),]

fit = lm(as.matrix(caray[,24:62]) ~ caray$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Alouatta_caraya.txt", col.names = NA)

# Alouatta fusca
fusca = msrs[which(msrs$SPECIES. == "fusca"),]

fit = lm(as.matrix(fusca[,24:62]) ~ fusca$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Alouatta_fusca.txt", col.names = NA)

# Alouatta palliata
palli = msrs[which(msrs$SPECIES. == "palliata"),]

fit = lm(as.matrix(palli[,24:62]) ~ palli$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Alouatta_palliata.txt", col.names = NA)

# Alouatta senicula
seni = msrs[which(msrs$SPECIES. == "senicula"),]

fit = lm(as.matrix(seni[,24:62]) ~ seni$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Alouatta_senicula.txt", col.names = NA)

# Alouatta villosa
villo = msrs[which(msrs$SPECIES. == "villosa"),]
# grupo irmao
palli = msrs[which(msrs$SPECIES. == "palliata"),]

villo = rbind(villo, palli)

fit = manova(as.matrix(villo[,24:62]) ~ villo$SEX4. + villo$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Alouatta_villosa.txt", col.names = NA)

# Aotus azarae
aza = msrs[which(msrs$SPECIES. == "azarae"),]
  
fit = lm(as.matrix(aza[,24:62]) ~ aza$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_azarae.txt", col.names = NA)

# Aotus brumbacki
brum = msrs[which(msrs$SPECIES. == "brumbacki"),]
# grupo irmao
voc = msrs[which(msrs$SPECIES. == "vociferans"),]
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]

brum = AppendMe(c("brum", "voc", "lem"))

fit = manova(as.matrix(brum[,24:62]) ~ brum$SEX4. + brum$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_brumbacki.txt", col.names = NA)

# Aotus infulatus
infu = msrs[which(msrs$SPECIES. == "infulatus"),]
  
fit = lm(as.matrix(infu[,24:62]) ~ infu$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_infulatus.txt", col.names = NA)

# Aotus lemurinus
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]

fit = manova(as.matrix(lem[,24:62]) ~ lem$SEX4. + lem$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_lemurinus.txt", col.names = NA)

# Aotus nancymai
nancy = msrs[which(msrs$SPECIES. == "nancymai"),]
# grupo irmao
aza = msrs[which(msrs$SPECIES. == "azarae"),]
nigri = msrs[which(msrs$SPECIES. == "nigriceps"),]

nancy = AppendMe(c("nancy", "aza", "nigri"))

fit = manova(as.matrix(nancy[,24:62]) ~ nancy$SEX4. + nancy$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_nancymai.txt", col.names = NA)

# Aotus nigriceps
nigri = msrs[which(msrs$SPECIES. == "nigriceps"),]

fit = lm(as.matrix(nigri[,24:62]) ~ nigri$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_nigriceps.txt", col.names = NA)

# Aotus trivirgatus
trivi = msrs[which(msrs$SPECIES. == "trivirgatus"),]
# grupo irmao
brum = msrs[which(msrs$SPECIES. == "brumbacki"),]
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]
vosc = msrs[which(msrs$SPECIES. == "vociferans"),]

trivi = AppendMe(c("trivi", "brum", "vosc", "lem"))

fit = manova(as.matrix(trivi[,24:62]) ~ trivi$SEX4. + trivi$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_trivirgatus.txt", col.names = NA)

# Aotus vociferans
vosc = msrs[which(msrs$SPECIES. == "vociferans"),]
# grupo irmao
brum = msrs[which(msrs$SPECIES. == "brumbacki"),]
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]

vosc = AppendMe(c("vosc", "brum", "lem"))

fit = manova(as.matrix(vosc[,24:62]) ~ vosc$SEX4. + vosc$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Aotus_vociferans.txt", col.names = NA)

# Ateles_belzebulth
belze = msrs[which(msrs$SPECIES. == "belzebulth"),]
# grupo irmao
geo = msrs[which(msrs$GENUS. == "Ateles"),]
geo = geo[which(geo$SPECIES. == "geoffroyi"),]

belze = rbind(belze, geo)

fit = manova(as.matrix(belze[,24:62]) ~ belze$SEX4. + belze$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Ateles_belzebulth.txt", col.names = NA)

# Ateles chamek
cham = msrs[which(msrs$SUB. == "chamek"),]
# grupo irmao
marg = msrs[which(msrs$SPECIES. == "marginatus"),]

cham = rbind(cham, marg)

fit = manova(as.matrix(cham[,24:62]) ~ cham$SEX4. + cham$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Ateles_chamek.txt", col.names = NA)

# Ateles geoffroyi
geo = msrs[which(msrs$GENUS. == "Ateles"),]
geo = geo[which(geo$SPECIES. == "geoffroyi"),]

fit = lm(as.matrix(geo[,24:62]) ~ geo$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Ateles_geoffroyi.txt", col.names = NA)

# Ateles marginatus
marg = msrs[which(msrs$SPECIES. == "marginatus"),]

fit = lm(as.matrix(marg[,24:62]) ~ marg$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Ateles_marginatus.txt", col.names = NA)

# Ateles paniscus
pan = msrs[which(msrs$SPECIES. == "paniscus"),]

fit = lm(as.matrix(pan[,24:62]) ~ pan$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Ateles_paniscus.txt", col.names = NA)

# Brachyteles arachnoides
brachy = msrs[which(msrs$SPECIES. == "arachnoides"),]

fit = lm(as.matrix(brachy[,24:62]) ~ brachy$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Brachyteles_arachnoides.txt", col.names = NA)

# Cacajao calvus
calvus = msrs[which(msrs$SPECIES. == "calvus"),]

fit = manova(as.matrix(calvus[,24:62]) ~ calvus$SEX4. + calvus$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cacajao_calvus.txt", col.names = NA)

# Cacajao melanocephalus
mela = msrs[which(msrs$SPECIES. == "melanocephalus"),]
# grupo irmao
calvus = msrs[which(msrs$SPECIES. == "calvus"),]

mela = rbind(mela, calvus)

fit = manova(as.matrix(mela[,24:62]) ~ mela$SEX4. + mela$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cacajao_melanocephalus.txt", col.names = NA)

# Callicebus brunneus
brun = msrs[which(msrs$SPECIES. == "brunneus"),]

fit = lm(as.matrix(brun[,24:62]) ~ brun$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callicebus_brunneus.txt", col.names = NA)

# Callicebus cupreus
cupre = msrs[which(msrs$SUB. == "cupreus"),]

fit = lm(as.matrix(cupre[,24:62]) ~ cupre$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callicebus_cupreus.txt", col.names = NA)

# Callicebus discolor
disco = msrs[which(msrs$SUB. == "discolor"),]

fit = lm(as.matrix(disco[,24:62]) ~ disco$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callicebus_discolor.txt", col.names = NA)

# Callicebus moloch
molo = msrs[which(msrs$SPECIES. == "moloch"),]

fit = lm(as.matrix(molo[,24:62]) ~ molo$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callicebus_moloch.txt", col.names = NA)

# Callicebus personatus
perso = msrs[which(msrs$SPECIES. == "personatus"),]

fit = manova(as.matrix(perso[,24:62]) ~ perso$SEX4. + perso$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callicebus_personatus.txt", col.names = NA)

# Callicebus torquatus
torq = msrs[which(msrs$SPECIES. == "torquatus"),]

fit = manova(as.matrix(torq[,24:62]) ~ torq$SEX4. + torq$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callicebus_torquatus.txt", col.names = NA)

# Callimico goeldii
goe = msrs[which(msrs$SUB. == "goeldii"),]

fit = lm(as.matrix(goe[,24:62]) ~ goe$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callimico_goeldii.txt", col.names = NA)

# Callithrix argentata
arge = msrs[which(msrs$SPECIES. == "argentata"),]

fit = lm(as.matrix(arge[,24:62]) ~ arge$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_argentata.txt", col.names = NA)

# Callithrix aurita
auri = msrs[which(msrs$SPECIES. == "aurita"),]
# grupo irmao
auri = msrs[which(msrs$GENUS. == "Callithrix"),]

fit = manova(as.matrix(auri[,24:62]) ~ auri$SEX4. + auri$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_aurita.txt", col.names = NA)

# Callithrix_emiliae
emi = msrs[which(msrs$SPECIES. == "emiliae"),]
# grupo irmao
sate = msrs[which(msrs$SPECIES. == "saterei"),]
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]
mau = msrs[which(msrs$SPECIES. == "mauesi"),]

emi = AppendMe(c("emi", "sate", "hume", "mau"))

fit = manova(as.matrix(emi[,24:62]) ~ emi$SEX4. + emi$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_emiliae.txt", col.names = NA)

# Callithrix geoffroyi
geo = msrs[which(msrs$GENUS. == "Callithrix"),]
geo = geo[which(geo$SPECIES. == "geoffroyi"),]

fit = lm(as.matrix(geo[,24:62]) ~ geo$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_geoffroyi.txt", col.names = NA)

# Callithrix humeralifera
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]

fit = lm(as.matrix(hume[,24:62]) ~ hume$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_humeralifera.txt", col.names = NA)

# Callithrix jacchus
jac = msrs[which(msrs$SPECIES. == "jacchus"),]

fit = manova(as.matrix(jac[,24:62]) ~ jac$SEX4. + jac$SUB. + 
               jac$SEX4.:jac$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_jacchus.txt", col.names = NA)

# Callithrix kuhlii
ku = msrs[which(msrs$SPECIES. == "kuhlii"),]

fit = lm(as.matrix(ku[,24:62]) ~ ku$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_kuhlii.txt", col.names = NA)

# Callithrix mauesi
mau = msrs[which(msrs$SPECIES. == "mauesi"),]
# grupo irmao
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]
sate = msrs[which(msrs$SPECIES. == "saterei"),]
emi = msrs[which(msrs$SPECIES. == "emiliae"),]

mau = AppendMe(c("mau", "hume", "sate", "emi"))

fit = manova(as.matrix(mau[,24:62]) ~ mau$SEX4. + mau$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_mauesi.txt", col.names = NA)

# Callithrix penicillata
peni = msrs[which(msrs$SPECIES. == "penicillata"),]

fit = lm(as.matrix(peni[,24:62]) ~ peni$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_penicillata.txt", col.names = NA)

# Callithrix saterei
sate = msrs[which(msrs$SPECIES. == "saterei"),]
# grupo irmao
mau = msrs[which(msrs$SPECIES. == "mauesi"),]
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]
emi = msrs[which(msrs$SPECIES. == "emiliae"),]

sate = AppendMe(c("sate", "mau", "hume", "emi"))

fit = manova(as.matrix(sate[,24:62]) ~ sate$SEX4. + sate$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Callithrix_saterei.txt", col.names = NA)

# Cebuella pygmaea
pyg = msrs[which(msrs$SPECIES. == "pygmaea"),]

fit = lm(as.matrix(pyg[,24:62]) ~ pyg$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebuella_pygmaea.txt", col.names = NA)

# Cebus albifrons
alb = msrs[which(msrs$SPECIES. == "albifrons"),]
# grupo irmao
cap = msrs[which(msrs$SPECIES. == "capucinus"),]
lib = msrs[which(msrs$SPECIES. == "libidinosus"),]
ape = msrs[which(msrs$SPECIES. == "apella"),]
rob = msrs[which(msrs$SPECIES. == "robustus"),]
xan1 = msrs[which(msrs$SPECIES. == "xanthosternos"),]
xan2 = msrs[which(msrs$SPECIES. == "xanthosternus"),]
xan = rbind(xan1,xan2)

alb = AppendMe(c("alb", "cap", "lib", "ape", "rob", "xan"))

fit = manova(as.matrix(alb[,24:62]) ~ alb$SEX4. + alb$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebus_albifrons.txt", col.names = NA)

# Cebus apella
ape = msrs[which(msrs$SPECIES. == "apella"),]

fit = lm(as.matrix(ape[,24:62]) ~ ape$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebus_apella.txt", col.names = NA)

# Cebus capucinus
cap = msrs[which(msrs$SPECIES. == "capucinus"),]
# grupo irmao
alb = msrs[which(msrs$SPECIES. == "albifrons"),]
lib = msrs[which(msrs$SPECIES. == "libidinosus"),]
ape = msrs[which(msrs$SPECIES. == "apella"),]
rob = msrs[which(msrs$SPECIES. == "robustus"),]
xan1 = msrs[which(msrs$SPECIES. == "xanthosternos"),]
xan2 = msrs[which(msrs$SPECIES. == "xanthosternus"),]
xan = rbind(xan1,xan2)

cap = AppendMe(c("alb", "cap", "lib", "ape", "rob", "xan"))

fit = manova(as.matrix(cap[,24:62]) ~ cap$SEX4. + cap$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebus_capucinus.txt", col.names = NA)

# Cebus libidinosus
lib = msrs[which(msrs$SPECIES. == "libidinosus"),]

fit = lm(as.matrix(lib[,24:62]) ~ lib$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebus_libidinosus.txt", col.names = NA)

# Cebus nigritus
nig = msrs[which(msrs$SPECIES. == "nigritus"),]

fit = manova(as.matrix(nig[,24:62]) ~ nig$SEX4. + nig$SUB. +
               nig$SEX4.:nig$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebus_nigritus.txt", col.names = NA)

# Cebus robustus
rob = msrs[which(msrs$SPECIES. == "robustus"),]

fit = lm(as.matrix(rob[,24:62]) ~ rob$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebus_robustus.txt", col.names = NA)

# Cebus xanthosternos
xan1 = msrs[which(msrs$SPECIES. == "xanthosternos"),]
xan2 = msrs[which(msrs$SPECIES. == "xanthosternus"),]
xan = rbind(xan1,xan2)
# grupo irmao
rob = msrs[which(msrs$SPECIES. == "robustus"),]
ape = msrs[which(msrs$SPECIES. == "apella"),]
lib = msrs[which(msrs$SPECIES. == "libidinosus"),]

xan = AppendMe(c("xan", "rob", "ape", "lib"))

fit = manova(as.matrix(xan[,24:62]) ~ xan$SEX4. + xan$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Cebus_xanthosternos.txt", col.names = NA)

# Chiropotes albinasus
albi = msrs[which(msrs$SPECIES. == "albinasus"),]

fit = lm(as.matrix(albi[,24:62]) ~ albi$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Chiropotes_albinasus.txt", col.names = NA)

# Chiropotes satanas
sata = msrs[which(msrs$SPECIES. == "satanas"),]

fit = manova(as.matrix(sata[,24:62]) ~ sata$SEX4. + sata$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Chiropotes_satanas.txt", col.names = NA)

# Lagothrix flavicauda
flavi = msrs[which(msrs$SPECIES. == "flavicauda"),]
# grupo irmao
tricha = msrs[which(msrs$SPECIES. == "lagothricha"),]

flavi = rbind(flavi, tricha)

fit = manova(as.matrix(flavi[,24:62]) ~ flavi$SEX4. + flavi$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Lagothrix_flavicauda.txt", col.names = NA)

# Lagothrix lagothricha
tricha = msrs[which(msrs$SPECIES. == "lagothricha"),]

fit = manova(as.matrix(tricha[,24:62]) ~ tricha$SEX4. + tricha$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Lagothrix_lagothricha.txt", col.names = NA)

# Leontopithecus chrysomelas
melas = msrs[which(msrs$SPECIES. == "chrysomelas"),]
# grupo irmao
rosa = msrs[which(msrs$SPECIES. == "rosalia"),]
pygus = msrs[which(msrs$SPECIES. == "chrysopygus"),]

melas = AppendMe(c("melas", "rosa", "pygus"))

fit = manova(as.matrix(melas[,24:62]) ~ melas$SEX4. + melas$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Leontopithecus_chrysomelas.txt", col.names = NA)

# Leontopithecus chrysopygus
pygus = msrs[which(msrs$SPECIES. == "chrysopygus"),]
# grupo irmao
rosa = msrs[which(msrs$SPECIES. == "rosalia"),]

pygus = rbind(pygus, rosa)

fit = manova(as.matrix(pygus[,24:62]) ~ pygus$SEX4. + pygus$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Leontopithecus_chrysopygus.txt", col.names = NA)

# Leontopithecus rosalia
rosa = msrs[which(msrs$SUB. == "rosalia"),]

fit = lm(as.matrix(rosa[,24:62]) ~ rosa$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Leontopithecus_rosalia.txt", col.names = NA)

# Pithecia irrorata
irro = msrs[which(msrs$SPECIES. == "irrorata"),]

fit = manova(as.matrix(irro[,24:62]) ~ irro$SEX4. + irro$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Pithecia_irrorata.txt", col.names = NA)

# Pithecia monacha
mona = msrs[which(msrs$SPECIES. == "monacha"),]

fit = lm(as.matrix(mona[,24:62]) ~ mona$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Pithecia_monacha.txt", col.names = NA)

# Pithecia pithecia
pite = msrs[which(msrs$SPECIES. == "pithecia"),]

fit = manova(as.matrix(pite[,24:62]) ~ pite$SEX4. + pite$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Pithecia_pithecia.txt", col.names = NA)

# Saguinus bicolor
bico = msrs[which(msrs$SPECIES. == "bicolor"),]
# grupo irmao
nig = msrs[which(msrs$SPECIES. == "midas"),]
nig = nig[which(nig$SUB. == "niger"),]

bico = rbind(bico, nig)

fit = manova(as.matrix(bico[,24:62]) ~ bico$SEX4. + bico$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_bicolor.txt", col.names = NA)

# Saguinus fuscicolis
fusci = msrs[which(msrs$SPECIES. == "fuscicollis"),]

fit = manova(as.matrix(fusci[,24:62]) ~ fusci$SEX4. + fusci$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_fuscicollis.txt", col.names = NA)

# Saguinus geoffroyi
geo = msrs[which(msrs$SPECIES. == "geoffroyi"),]
geo = geo[which(geo$GENUS. == "Saguinus"),]

fit = lm(as.matrix(geo[,24:62]) ~ geo$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_geoffroyi.txt", col.names = NA)

# Saguinus imperator
impe = msrs[which(msrs$SPECIES. == "imperator"),]

fit = lm(as.matrix(impe[,24:62]) ~ impe$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_imperator.txt", col.names = NA)

# Saguinus labiatus
labi = msrs[which(msrs$SPECIES. == "labiatus"),]

fit = lm(as.matrix(labi[,24:62]) ~ labi$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_labiatus.txt", col.names = NA)

# Saguinus leucopus
leuco = msrs[which(msrs$SPECIES. == "leucopus"),]

fit = lm(as.matrix(leuco[,24:62]) ~ leuco$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_leucopus.txt", col.names = NA)

# Saguinus midas
midas = msrs[which(msrs$SPECIES. == "midas"),]

fit = manova(as.matrix(midas[,24:62]) ~ midas$SEX4. + midas$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_midas.txt", col.names = NA)

# Saguinus mystax
mis = msrs[which(msrs$SPECIES. == "mystax"),]

fit = manova(as.matrix(mis[,24:62]) ~ mis$SEX4. + mis$SUB. +
               mis$SEX4.:mis$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_mystax.txt", col.names = NA)

# Saguinus nigricollis
nigri = msrs[which(msrs$SPECIES. == "nigricollis"),]

fit = manova(as.matrix(nigri[,24:62]) ~ nigri$SEX4. + nigri$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_nigricollis.txt", col.names = NA)

# Saguinus oedipus
oedi = msrs[which(msrs$SPECIES. == "oedipus"),]

fit = lm(as.matrix(oedi[,24:62]) ~ oedi$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saguinus_oedipus.txt", col.names = NA)

# Saimiri boliviensis
boli = msrs[which(msrs$SPECIES. == "boliviensis"),]
# hibridos
x = msrs[which(msrs$SPECIES. == "boliviensis x macrodon"),]

boli = rbind(boli, x)

fit = lm(as.matrix(boli[,24:62]) ~ boli$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saimiri_boliviensis.txt", col.names = NA)

# Saimiri cassiquiarensis
cassi = msrs[which(msrs$SPECIES. == "cassiquiarensis"),]
# hibridos
x = msrs[which(msrs$SPECIES. == "cassiquiarensis x sciureus"),]

cassi = rbind(cassi, x)

fit = lm(as.matrix(cassi[,24:62]) ~ cassi$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saimiri_cassiquiarensis.txt", col.names = NA)

# Saimiri oerstedi
oer = msrs[which(msrs$SPECIES. == "oerstedi"),]

fit = lm(as.matrix(oer[,24:62]) ~ oer$SEX4.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saimiri_oerstedi.txt", col.names = NA)

# Saimiri sciureus
sci = msrs[which(msrs$SPECIES. == "sciureus"),]
# hibridos
x1 = msrs[which(msrs$SPECIES. == "sciureus x cassiquiarensis"),]
y1 = msrs[which(msrs$SPECIES. == "sciureus x ustus"),]
z1 = msrs[which(msrs$SPECIES. == "sciurues x cassiquiarensis"),]

sci = AppendMe(c("sci", "x1", "y1", "z1"))

fit = manova(as.matrix(sci[,24:62]) ~ sci$SEX4. + sci$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saimiri_sciureus.txt", col.names = NA)

# Saimiri ustus
ustu = msrs[which(msrs$SPECIES. == "ustus"),]
# hibridos
x1 = msrs[which(msrs$SPECIES. == "ustus x macrodon"),]
y1 = msrs[which(msrs$SPECIES. == "ustus x sciureus"),]

ustu = AppendMe(c("ustu", "x1", "y1"))

fit = manova(as.matrix(ustu[,24:62]) ~ ustu$SEX4. + ustu$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_nolog_vcv/Saimiri_ustus.txt", col.names = NA)


