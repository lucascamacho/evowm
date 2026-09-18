# Norms of the vector dimorphism
# set wd

# vector correlation vector
prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# load dimorphism data
averages = readRDS("~/Dropbox/Novo_Output/Averages_PCS_Ancestral_Species.RData")

# get dimorphism data and apply Norma function
nolog_dimorphism = averages$NoLog_ByTrait_Ancestral_Dimorphism
log_dimorphism = averages$Log_ByTrait_Ancestral_Dimorphism
  
nolog_norms = apply(nolog_dimorphism, 1, norma)
log_norms = apply(log_dimorphism, 1, norma)

# save results
normas = data.frame(nolog_norms, log_norms)
colnames(normas) = c("NoLog_Norms", "Log_Norms")

saveRDS(normas, file = "~/Dropbox/Novo_Output/Ancestral_Dimorphism_Norms.RData")
