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

fit = manova(as.matrix(log(belze[,24:62])) ~ belze$SEX. + belze$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Alouatta_belzebul.txt", col.names = NA)

# Alouatta caraya
caray = msrs[which(msrs$SPECIES. == "caraya"),]

fit = lm(as.matrix(log(caray[,24:62])) ~ caray$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Alouatta_caraya.txt", col.names = NA)

# Alouatta fusca
fusca = msrs[which(msrs$SPECIES. == "fusca"),]

fit = lm(as.matrix(log(fusca[,24:62])) ~ fusca$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Alouatta_fusca.txt", col.names = NA)

# Alouatta palliata
palli = msrs[which(msrs$SPECIES. == "palliata"),]

fit = lm(as.matrix(log(palli[,24:62])) ~ palli$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Alouatta_palliata.txt", col.names = NA)

# Alouatta senicula
seni = msrs[which(msrs$SPECIES. == "senicula"),]

fit = lm(as.matrix(log(seni[,24:62])) ~ seni$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Alouatta_senicula.txt", col.names = NA)

# Alouatta villosa
villo = msrs[which(msrs$SPECIES. == "villosa"),]
# grupo irmao
palli = msrs[which(msrs$SPECIES. == "palliata"),]

villo = rbind(villo, palli)

fit = manova(as.matrix(log(villo[,24:62])) ~ villo$SEX. + villo$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Alouatta_villosa.txt", col.names = NA)

# Aotus azarae
aza = msrs[which(msrs$SPECIES. == "azarae"),]

fit = lm(as.matrix(log(aza[,24:62])) ~ aza$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_azarae.txt", col.names = NA)

# Aotus brumbacki
brum = msrs[which(msrs$SPECIES. == "brumbacki"),]
# grupo irmao
voc = msrs[which(msrs$SPECIES. == "vociferans"),]
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]

brum = AppendMe(c("brum", "voc", "lem"))

fit = manova(as.matrix(log(brum[,24:62])) ~ brum$SEX. + brum$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_brumbacki.txt", col.names = NA)

# Aotus infulatus
infu = msrs[which(msrs$SPECIES. == "infulatus"),]

fit = lm(as.matrix(log(infu[,24:62])) ~ infu$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_infulatus.txt", col.names = NA)

# Aotus lemurinus
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]

fit = manova(as.matrix(log(lem[,24:62])) ~ lem$SEX. + lem$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_lemurinus.txt", col.names = NA)

# Aotus nancymai
nancy = msrs[which(msrs$SPECIES. == "nancymai"),]
# grupo irmao
aza = msrs[which(msrs$SPECIES. == "azarae"),]
nigri = msrs[which(msrs$SPECIES. == "nigriceps"),]

nancy = AppendMe(c("nancy", "aza", "nigri"))

fit = manova(as.matrix(log(nancy[,24:62])) ~ nancy$SEX. + nancy$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_nancymai.txt", col.names = NA)

# Aotus nigriceps
nigri = msrs[which(msrs$SPECIES. == "nigriceps"),]

fit = lm(as.matrix(log(nigri[,24:62])) ~ nigri$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_nigriceps.txt", col.names = NA)

# Aotus trivirgatus
trivi = msrs[which(msrs$SPECIES. == "trivirgatus"),]
# grupo irmao
brum = msrs[which(msrs$SPECIES. == "brumbacki"),]
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]
vosc = msrs[which(msrs$SPECIES. == "vociferans"),]

trivi = AppendMe(c("trivi", "brum", "vosc", "lem"))

fit = manova(as.matrix(log(trivi[,24:62])) ~ trivi$SEX. + trivi$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_trivirgatus.txt", col.names = NA)

# Aotus vociferans
vosc = msrs[which(msrs$SPECIES. == "vociferans"),]
# grupo irmao
brum = msrs[which(msrs$SPECIES. == "brumbacki"),]
lem = msrs[which(msrs$SPECIES. == "lemurinus"),]

vosc = AppendMe(c("vosc", "brum", "lem"))

fit = manova(as.matrix(log(vosc[,24:62])) ~ vosc$SEX. + vosc$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Aotus_vociferans.txt", col.names = NA)

# Ateles_belzebulth
belze = msrs[which(msrs$SPECIES. == "belzebulth"),]
# grupo irmao
geo = msrs[which(msrs$GENUS. == "Ateles"),]
geo = geo[which(geo$SPECIES. == "geoffroyi"),]

belze = rbind(belze, geo)

fit = manova(as.matrix(log(belze[,24:62])) ~ belze$SEX. + belze$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Ateles_belzebulth.txt", col.names = NA)

# Ateles chamek
cham = msrs[which(msrs$SUB. == "chamek"),]
# grupo irmao
marg = msrs[which(msrs$SPECIES. == "marginatus"),]

cham = rbind(cham, marg)

fit = manova(as.matrix(log(cham[,24:62])) ~ cham$SEX. + cham$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Ateles_chamek.txt", col.names = NA)

# Ateles geoffroyi
geo = msrs[which(msrs$GENUS. == "Ateles"),]
geo = geo[which(geo$SPECIES. == "geoffroyi"),]

fit = lm(as.matrix(log(geo[,24:62])) ~ geo$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Ateles_geoffroyi.txt", col.names = NA)

# Ateles marginatus
marg = msrs[which(msrs$SPECIES. == "marginatus"),]

fit = lm(as.matrix(log(marg[,24:62])) ~ marg$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Ateles_marginatus.txt", col.names = NA)

# Ateles paniscus
pan = msrs[which(msrs$SPECIES. == "paniscus"),]

fit = lm(as.matrix(log(pan[,24:62])) ~ pan$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Ateles_paniscus.txt", col.names = NA)

# Brachyteles arachnoides
brachy = msrs[which(msrs$SPECIES. == "arachnoides"),]

fit = lm(as.matrix(log(brachy[,24:62])) ~ brachy$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Brachyteles_arachnoides.txt", col.names = NA)

# Cacajao calvus
calvus = msrs[which(msrs$SPECIES. == "calvus"),]

fit = manova(as.matrix(log(calvus[,24:62])) ~ calvus$SEX. + calvus$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cacajao_calvus.txt", col.names = NA)

# Cacajao melanocephalus
mela = msrs[which(msrs$SPECIES. == "melanocephalus"),]
# grupo irmao
calvus = msrs[which(msrs$SPECIES. == "calvus"),]

mela = rbind(mela, calvus)

fit = manova(as.matrix(log(mela[,24:62])) ~ mela$SEX. + mela$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cacajao_melanocephalus.txt", col.names = NA)

# Callicebus brunneus
brun = msrs[which(msrs$SPECIES. == "brunneus"),]

fit = lm(as.matrix(log(brun[,24:62])) ~ brun$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callicebus_brunneus.txt", col.names = NA)

# Callicebus cupreus
cupre = msrs[which(msrs$SUB. == "cupreus"),]

fit = lm(as.matrix(log(cupre[,24:62])) ~ cupre$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callicebus_cupreus.txt", col.names = NA)

# Callicebus discolor
disco = msrs[which(msrs$SUB. == "discolor"),]

fit = lm(as.matrix(log(disco[,24:62])) ~ disco$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callicebus_discolor.txt", col.names = NA)

# Callicebus moloch
molo = msrs[which(msrs$SPECIES. == "moloch"),]

fit = lm(as.matrix(log(molo[,24:62])) ~ molo$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callicebus_moloch.txt", col.names = NA)

# Callicebus personatus
perso = msrs[which(msrs$SPECIES. == "personatus"),]

fit = manova(as.matrix(log(perso[,24:62])) ~ perso$SEX. + perso$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callicebus_personatus.txt", col.names = NA)

# Callicebus torquatus
torq = msrs[which(msrs$SPECIES. == "torquatus"),]

fit = manova(as.matrix(log(torq[,24:62])) ~ torq$SEX. + torq$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callicebus_torquatus.txt", col.names = NA)

# Callimico goeldii
goe = msrs[which(msrs$SUB. == "goeldii"),]

fit = lm(as.matrix(log(goe[,24:62])) ~ goe$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callimico_goeldii.txt", col.names = NA)

# Callithrix argentata
arge = msrs[which(msrs$SPECIES. == "argentata"),]

fit = lm(as.matrix(log(arge[,24:62])) ~ arge$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_argentata.txt", col.names = NA)

# Callithrix aurita
auri = msrs[which(msrs$SPECIES. == "aurita"),]
# grupo irmao
auri = msrs[which(msrs$GENUS. == "Callithrix"),]

fit = manova(as.matrix(log(auri[,24:62])) ~ auri$SEX. + auri$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_aurita.txt", col.names = NA)

# Callithrix_emiliae
emi = msrs[which(msrs$SPECIES. == "emiliae"),]
# grupo irmao
sate = msrs[which(msrs$SPECIES. == "saterei"),]
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]
mau = msrs[which(msrs$SPECIES. == "mauesi"),]

emi = AppendMe(c("emi", "sate", "hume", "mau"))

fit = manova(as.matrix(log(emi[,24:62])) ~ emi$SEX. + emi$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_emiliae.txt", col.names = NA)

# Callithrix geoffroyi
geo = msrs[which(msrs$GENUS. == "Callithrix"),]
geo = geo[which(geo$SPECIES. == "geoffroyi"),]

fit = lm(as.matrix(log(geo[,24:62])) ~ geo$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_geoffroyi.txt", col.names = NA)

# Callithrix humeralifera
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]

fit = lm(as.matrix(log(hume[,24:62])) ~ hume$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_humeralifera.txt", col.names = NA)

# Callithrix jacchus
jac = msrs[which(msrs$SPECIES. == "jacchus"),]

fit = manova(as.matrix(log(jac[,24:62])) ~ jac$SEX. + jac$SUB. + 
               jac$SEX.:jac$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_jacchus.txt", col.names = NA)

# Callithrix kuhlii
ku = msrs[which(msrs$SPECIES. == "kuhlii"),]

fit = lm(as.matrix(log(ku[,24:62])) ~ ku$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_kuhlii.txt", col.names = NA)

# Callithrix mauesi
mau = msrs[which(msrs$SPECIES. == "mauesi"),]
# grupo irmao
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]
sate = msrs[which(msrs$SPECIES. == "saterei"),]
emi = msrs[which(msrs$SPECIES. == "emiliae"),]

mau = AppendMe(c("mau", "hume", "sate", "emi"))

fit = manova(as.matrix(log(mau[,24:62])) ~ mau$SEX4. + mau$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_mauesi.txt", col.names = NA)

# Callithrix penicillata
peni = msrs[which(msrs$SPECIES. == "penicillata"),]

fit = lm(as.matrix(log(peni[,24:62])) ~ peni$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_penicillata.txt", col.names = NA)

# Callithrix saterei
sate = msrs[which(msrs$SPECIES. == "saterei"),]
# grupo irmao
mau = msrs[which(msrs$SPECIES. == "mauesi"),]
hume = msrs[which(msrs$SPECIES. == "humeralifera"),]
emi = msrs[which(msrs$SPECIES. == "emiliae"),]

sate = AppendMe(c("sate", "mau", "hume", "emi"))

fit = manova(as.matrix(log(sate[,24:62])) ~ sate$SEX4. + sate$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Callithrix_saterei.txt", col.names = NA)

# Cebuella pygmaea
pyg = msrs[which(msrs$SPECIES. == "pygmaea"),]

fit = lm(as.matrix(log(pyg[,24:62])) ~ pyg$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebuella_pygmaea.txt", col.names = NA)

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

fit = manova(as.matrix(log(alb[,24:62])) ~ alb$SEX4. + alb$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebus_albifrons.txt", col.names = NA)

# Cebus apella
ape = msrs[which(msrs$SPECIES. == "apella"),]

fit = lm(as.matrix(log(ape[,24:62])) ~ ape$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebus_apella.txt", col.names = NA)

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

fit = manova(as.matrix(log(cap[,24:62])) ~ cap$SEX4. + cap$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebus_capucinus.txt", col.names = NA)

# Cebus libidinosus
lib = msrs[which(msrs$SPECIES. == "libidinosus"),]

fit = lm(as.matrix(log(lib[,24:62])) ~ lib$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebus_libidinosus.txt", col.names = NA)

# Cebus nigritus
nig = msrs[which(msrs$SPECIES. == "nigritus"),]

fit = manova(as.matrix(log(nig[,24:62])) ~ nig$SEX. + nig$SUB. +
               nig$SEX.:nig$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebus_nigritus.txt", col.names = NA)

# Cebus robustus
rob = msrs[which(msrs$SPECIES. == "robustus"),]

fit = lm(as.matrix(log(rob[,24:62])) ~ rob$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebus_robustus.txt", col.names = NA)

# Cebus xanthosternos
xan1 = msrs[which(msrs$SPECIES. == "xanthosternos"),]
xan2 = msrs[which(msrs$SPECIES. == "xanthosternus"),]
xan = rbind(xan1,xan2)
# grupo irmao
rob = msrs[which(msrs$SPECIES. == "robustus"),]
ape = msrs[which(msrs$SPECIES. == "apella"),]
lib = msrs[which(msrs$SPECIES. == "libidinosus"),]

xan = AppendMe(c("xan", "rob", "ape", "lib"))

fit = manova(as.matrix(log(xan[,24:62])) ~ xan$SEX. + xan$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Cebus_xanthosternos.txt", col.names = NA)

# Chiropotes albinasus
albi = msrs[which(msrs$SPECIES. == "albinasus"),]

fit = lm(as.matrix(log(albi[,24:62])) ~ albi$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Chiropotes_albinasus.txt", col.names = NA)

# Chiropotes satanas
sata = msrs[which(msrs$SPECIES. == "satanas"),]

fit = manova(as.matrix(log(sata[,24:62])) ~ sata$SEX. + sata$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Chiropotes_satanas.txt", col.names = NA)

# Lagothrix flavicauda
flavi = msrs[which(msrs$SPECIES. == "flavicauda"),]
# grupo irmao
tricha = msrs[which(msrs$SPECIES. == "lagothricha"),]

flavi = rbind(flavi, tricha)

fit = manova(as.matrix(log(flavi[,24:62])) ~ flavi$SEX. + flavi$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Lagothrix_flavicauda.txt", col.names = NA)

# Lagothrix lagothricha
tricha = msrs[which(msrs$SPECIES. == "lagothricha"),]

fit = manova(as.matrix(log(tricha[,24:62])) ~ tricha$SEX. + tricha$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Lagothrix_lagothricha.txt", col.names = NA)

# Leontopithecus chrysomelas
melas = msrs[which(msrs$SPECIES. == "chrysomelas"),]
# grupo irmao
rosa = msrs[which(msrs$SPECIES. == "rosalia"),]
pygus = msrs[which(msrs$SPECIES. == "chrysopygus"),]

melas = AppendMe(c("melas", "rosa", "pygus"))

fit = manova(as.matrix(log(melas[,24:62])) ~ melas$SEX. + melas$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Leontopithecus_chrysomelas.txt", col.names = NA)

# Leontopithecus chrysopygus
pygus = msrs[which(msrs$SPECIES. == "chrysopygus"),]
# grupo irmao
rosa = msrs[which(msrs$SPECIES. == "rosalia"),]

pygus = rbind(pygus, rosa)

fit = manova(as.matrix(log(pygus[,24:62])) ~ pygus$SEX. + pygus$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Leontopithecus_chrysopygus.txt", col.names = NA)

# Leontopithecus rosalia
rosa = msrs[which(msrs$SUB. == "rosalia"),]

fit = lm(as.matrix(log(rosa[,24:62])) ~ rosa$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Leontopithecus_rosalia.txt", col.names = NA)

# Pithecia irrorata
irro = msrs[which(msrs$SPECIES. == "irrorata"),]

fit = manova(as.matrix(log(irro[,24:62])) ~ irro$SEX. + irro$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Pithecia_irrorata.txt", col.names = NA)

# Pithecia monacha
mona = msrs[which(msrs$SPECIES. == "monacha"),]

fit = lm(as.matrix(log(mona[,24:62])) ~ mona$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Pithecia_monacha.txt", col.names = NA)

# Pithecia pithecia
pite = msrs[which(msrs$SPECIES. == "pithecia"),]

fit = manova(as.matrix(log(pite[,24:62])) ~ pite$SEX. + pite$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Pithecia_pithecia.txt", col.names = NA)

# Saguinus bicolor
bico = msrs[which(msrs$SPECIES. == "bicolor"),]
# grupo irmao
nig = msrs[which(msrs$SPECIES. == "midas"),]
nig = nig[which(nig$SUB. == "niger"),]

bico = rbind(bico, nig)

fit = manova(as.matrix(log(bico[,24:62])) ~ bico$SEX. + bico$SPECIES.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_bicolor.txt", col.names = NA)

# Saguinus fuscicolis
fusci = msrs[which(msrs$SPECIES. == "fuscicollis"),]

fit = manova(as.matrix(log(fusci[,24:62])) ~ fusci$SEX. + fusci$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_fuscicollis.txt", col.names = NA)

# Saguinus geoffroyi
geo = msrs[which(msrs$SPECIES. == "geoffroyi"),]
geo = geo[which(geo$GENUS. == "Saguinus"),]

fit = lm(as.matrix(log(geo[,24:62])) ~ geo$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_geoffroyi.txt", col.names = NA)

# Saguinus imperator
impe = msrs[which(msrs$SPECIES. == "imperator"),]

fit = lm(as.matrix(log(impe[,24:62])) ~ impe$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_imperator.txt", col.names = NA)

# Saguinus labiatus
labi = msrs[which(msrs$SPECIES. == "labiatus"),]

fit = lm(as.matrix(log(labi[,24:62])) ~ labi$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_labiatus.txt", col.names = NA)

# Saguinus leucopus
leuco = msrs[which(msrs$SPECIES. == "leucopus"),]

fit = lm(as.matrix(log(leuco[,24:62])) ~ leuco$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_leucopus.txt", col.names = NA)

# Saguinus midas
midas = msrs[which(msrs$SPECIES. == "midas"),]

fit = manova(as.matrix(log(midas[,24:62])) ~ midas$SEX. + midas$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_midas.txt", col.names = NA)

# Saguinus mystax
mis = msrs[which(msrs$SPECIES. == "mystax"),]

fit = manova(as.matrix(log(mis[,24:62])) ~ mis$SEX. + mis$SUB. +
               mis$SEX.:mis$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_mystax.txt", col.names = NA)

# Saguinus nigricollis
nigri = msrs[which(msrs$SPECIES. == "nigricollis"),]

fit = manova(as.matrix(log(nigri[,24:62])) ~ nigri$SEX. + nigri$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_nigricollis.txt", col.names = NA)

# Saguinus oedipus
oedi = msrs[which(msrs$SPECIES. == "oedipus"),]

fit = lm(as.matrix(log(oedi[,24:62])) ~ oedi$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saguinus_oedipus.txt", col.names = NA)

# Saimiri boliviensis
boli = msrs[which(msrs$SPECIES. == "boliviensis"),]
# hibridos
x = msrs[which(msrs$SPECIES. == "boliviensis x macrodon"),]

boli = rbind(boli, x)

fit = lm(as.matrix(log(boli[,24:62])) ~ boli$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saimiri_boliviensis.txt", col.names = NA)

# Saimiri cassiquiarensis
cassi = msrs[which(msrs$SPECIES. == "cassiquiarensis"),]
# hibridos
x = msrs[which(msrs$SPECIES. == "cassiquiarensis x sciureus"),]

cassi = rbind(cassi, x)

fit = lm(as.matrix(log(cassi[,24:62])) ~ cassi$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saimiri_cassiquiarensis.txt", col.names = NA)

# Saimiri oerstedi
oer = msrs[which(msrs$SPECIES. == "oerstedi"),]

fit = lm(as.matrix(log(oer[,24:62])) ~ oer$SEX.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saimiri_oerstedi.txt", col.names = NA)

# Saimiri sciureus
sci = msrs[which(msrs$SPECIES. == "sciureus"),]
# hibridos
x1 = msrs[which(msrs$SPECIES. == "sciureus x cassiquiarensis"),]
y1 = msrs[which(msrs$SPECIES. == "sciureus x ustus"),]
z1 = msrs[which(msrs$SPECIES. == "sciurues x cassiquiarensis"),]

sci = AppendMe(c("sci", "x1", "y1", "z1"))

fit = manova(as.matrix(log(sci[,24:62])) ~ sci$SEX. + sci$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saimiri_sciureus.txt", col.names = NA)

# Saimiri ustus
ustu = msrs[which(msrs$SPECIES. == "ustus"),]
# hibridos
x1 = msrs[which(msrs$SPECIES. == "ustus x macrodon"),]
y1 = msrs[which(msrs$SPECIES. == "ustus x sciureus"),]

ustu = AppendMe(c("ustu", "x1", "y1"))

fit = manova(as.matrix(log(ustu[,24:62])) ~ ustu$SEX. + ustu$SUB.)

summary(fit, test = "Wilks")
cov.matrix = CalculateMatrix(fit)

write.table(cov.matrix, file = "~/Dropbox/Doc/Output/p_log_vcv/Saimiri_ustus.txt", col.names = NA)


