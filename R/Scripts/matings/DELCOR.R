############################################################
# Lagged evolutionary correlation (TNT-like) in R
############################################################

library(ape)
library(phytools)
library(dplyr)
library(ggplot2)

#############################
# 1. INPUTS
#############################
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/")

library(ape)
library(dplyr)
library(phytools)

# carregando dados
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
ancestrals = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/ancestrals_averages_PCS_autovalues_primates.RDS")
#nicho = readxl::read_xls("~/Dropbox/Doc/Data/nicho.xls")

# pegar especies comuns entre nicho e o resto dos dados
sp_comuns <- Reduce(intersect, list(
  #nicho$genus,
  align$matings.genus,
  evolvas$genus,
  mds$genus
))

#nicho_f   <- nicho[nicho$genus %in% sp_comuns, ]
align_f   <- align[align$matings.genus %in% sp_comuns, ]
evolvas_f <- evolvas[evolvas$genus %in% sp_comuns, ]
mds_f     <- mds[mds$genus %in% sp_comuns, ]

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

# read and plot phylo tree
commom <- intersect(names(vcv), mds_f$genus)
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = ape::read.nexus(filename)
genus <- sub("_.*", "", tree$tip.label)
keep <- genus %in% commom
tree2 <- drop.tip(tree, tree$tip.label[!keep])

tree <- tree2
# extrair generos
genera <- sapply(strsplit(tree$tip.label, "_"), `[`, 1)
unique_genera <- unique(genera)
representatives <- c()
for(g in unique_genera){
  
  spp <- tree$tip.label[genera == g]
  
  if(length(spp) == 1){
    
    representatives[g] <- spp
    
  } else {
    
    # pegar comprimento do ramo terminal
    edges <- match(spp, tree$tip.label)
    terminal_lengths <- tree$edge.length[match(edges, tree$edge[,2])]
    
    # escolher o menor (mais recente)
    representatives[g] <- spp[which.max(terminal_lengths)]
  }
}

# remover todas as outras espécies
tree_genus <- drop.tip(tree, setdiff(tree$tip.label, representatives))

# renomear para os gêneros
tree_genus$tip.label <- names(representatives)

# #
# media_genero <- nicho %>%
#   dplyr::group_by(.data$genus) %>%
#   dplyr::summarise(
#     climacv_mean = mean(.data$climacv, na.rm = TRUE),
#     CVn_mean = mean(.data$CVn, na.rm = TRUE)
#   )

#nicho_f <- media_genero[match(tree_genus$tip.label, media_genero$genus), ]
evolvas_f <- evolvas_f[match(tree_genus$tip.label, evolvas_f$genus), ]
align_f <- align_f[match(tree_genus$tip.label, align_f$matings.genus), ]
mds_f <- mds_f[match(tree_genus$tip.label, mds_f$genus), ]

################################################################################
# HAPLORRHINI
df <- data.frame(
  genus = mds_f$genus,
  parvorder = mds_f$dados.PARVORDER,
  social = mds_f$dados.SOCIAL_ORGANIZATION,
  mating = mds_f$dados.MATING_SYSTEM,
  prop = mds_f$dados.PROP_MALES_FEMALES,
  dominance = mds_f$dados.DOMINANCE,
  aggression = mds_f$dados.AGGRESSION,
  dimorfism = evolvas_f$Dimorfism,
  integration = evolvas_f$Integration,
  size = evolvas_f$Size,
  nmds1 = mds_f$vegan..scores.fit..sites...1.,
  nmds2 = mds_f$vegan..scores.fit..sites...2.,
  normas = align_f$normas
  #clima = nicho_f$climacv_mean,
  #cvn = nicho_f$CVn_mean
)

tree <- tree_genus
# X <- evolvas_f$Integration
# Y <- align_f$normas
# names(X) <- evolvas_f$genus
# names(Y) <- align_f$matings.genus

# X <- evolvas_f$Size
# Y <- align_f$normas
# names(X) <- evolvas_f$genus
# names(Y) <- align_f$matings.genus

# X <- df$nmds1
# Y <- align_f$normas
# names(X) <- df$genus
# names(Y) <- align_f$matings.genus

X <- df$nmds2
Y <- align_f$normas
names(X) <- df$genus
names(Y) <- align_f$matings.genus
################################################################################
# CATARRHINI
align_c = align_f[1:16,]
evolvas_c = evolvas_f[1:16,]
mds_c = mds_f[1:16,]
tree_c = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_c$matings.genus))
nicho_c = nicho_f[1:16,]

df <- data.frame(
  genus = mds_c$genus,
  parvorder = mds_c$dados.PARVORDER,
  social = mds_c$dados.SOCIAL_ORGANIZATION,
  mating = mds_c$dados.MATING_SYSTEM,
  prop = mds_c$dados.PROP_MALES_FEMALES,
  dominance = mds_c$dados.DOMINANCE,
  aggression = mds_c$dados.AGGRESSION,
  dimorfism = evolvas_c$Dimorfism,
  integration = evolvas_c$Integration,
  size = evolvas_c$Size,
  nmds1 = mds_c$vegan..scores.fit..sites...1.,
  nmds2 = mds_c$vegan..scores.fit..sites...2.,
  normas = align_c$normas,
  clima = nicho_c$climacv_mean,
  cvn = nicho_c$CVn_mean
)

df <- df[match(tree_c$tip.label, df$genus), ]

tree <- tree_c
#X <- align_c$normas
#Y <- evolvas_c$Integration
#names(X) <- align_c$matings.genus
#names(Y) <- evolvas_c$genus
X <- evolvas_c$Integration
Y <- align_c$normas
names(X) <- evolvas_c$genus
names(Y) <- align_c$matings.genus

################################################################################
# PLATYRRHINI
align_p = align_f[17:28,]
evolvas_p = evolvas_f[17:28,]
mds_p = mds_f[17:28,]
tree_p = drop.tip(tree_genus, setdiff(tree_genus$tip.label, align_p$matings.genus))
nicho_p = nicho_f[17:28,]

#
df <- data.frame(
  genus = mds_p$genus,
  parvorder = mds_p$dados.PARVORDER,
  social = mds_p$dados.SOCIAL_ORGANIZATION,
  mating = mds_p$dados.MATING_SYSTEM,
  prop = mds_p$dados.PROP_MALES_FEMALES,
  dominance = mds_p$dados.DOMINANCE,
  aggression = mds_p$dados.AGGRESSION,
  dimorfism = evolvas_p$Dimorfism,
  integration = evolvas_p$Integration,
  size = evolvas_p$Size,
  nmds1 = mds_p$vegan..scores.fit..sites...1.,
  nmds2 = mds_p$vegan..scores.fit..sites...2.,
  normas = align_p$normas,
  clima = nicho_p$climacv_mean,
  cvn = nicho_p$CVn_mean
)

df <- df[match(tree_p$tip.label, df$genus), ]

tree <- tree_p
#X <- align_p$normas
#Y <- evolvas_p$Integration
#names(X) <- align_p$matings.genus
#names(Y) <- evolvas_p$genus
X <- evolvas_p$Integration
Y <- align_p$normas
names(X) <- evolvas_p$genus
names(Y) <- align_p$matings.genus

###############################################################################
## DELCOR
###############################################################################
anc_X <- fastAnc(tree, X)
anc_Y <- fastAnc(tree, Y)

# Combine tip + internal node states
states_X <- c(X, anc_X)
states_Y <- c(Y, anc_Y)

states_X <- states_X[tree$tip.label]
states_Y <- states_Y[tree$tip.label]

Ntip <- length(tree$tip.label)
Nnode <- tree$Nnode

#############################
# 3. COMPUTE BRANCH INCREMENTS
#############################
get_increments <- function(tree, states){
  parent <- as.character(tree$edge[,1])
  child  <- as.character(tree$edge[,2])
  
  states[child] - states[parent]
}

# X
states_vec <- numeric(Ntip + Nnode)
names(states_vec) <- as.character(1:(Ntip + Nnode))
states_vec[as.character(1:Ntip)] <- states_X[tree$tip.label]
states_vec[as.character((Ntip + 1):(Ntip + Nnode))] <- anc_X
incX <- get_increments(tree, states_vec)

# Y
states_vec <- numeric(Ntip + Nnode)
names(states_vec) <- as.character(1:(Ntip + Nnode))
states_vec[as.character(1:Ntip)] <- states_Y[tree$tip.label]
states_vec[as.character((Ntip + 1):(Ntip + Nnode))] <- anc_Y
incY <- get_increments(tree, states_vec)

#############################
# 4. FIND LAGGED PAIRS (TNT-like logic)
#############################

find_pairs <- function(tree, incx, incy,
                       radius = 4,
                       minx = 1e-4,
                       miny = 1e-4,
                       delfac = 1){
  
  pairs <- data.frame(x = numeric(),
                      y = numeric(),
                      w = numeric())
  
  for(i in 1:length(incx)){
    
    xchange <- incx[i]
    
    if(!is.na(xchange) && abs(xchange) >= minx){
      
      node <- tree$edge[i,2]
      dist <- 0
      
      repeat{
        
        parent_edge <- which(tree$edge[,2] == node)
        
        if(length(parent_edge) == 0) break
        
        parent <- tree$edge[parent_edge,1]
        dist <- dist + 1
        
        ychange <- incy[parent_edge]
        
        if(!is.na(ychange) && (abs(ychange) >= miny || dist > radius)){
          
          w <- 1 / (1 + dist * delfac)
          
          pairs <- rbind(pairs,
                         data.frame(x = xchange,
                                    y = ychange,
                                    w = w))
          
          break
        }
        
        node <- parent
      }
      
    }
    
  }
  
  return(pairs)
}

#############################
# 5. EXTRACT PAIRS
#############################

pairs <- find_pairs(tree, incX, incY)

# Remove any NA pairs just in case
pairs <- na.omit(pairs)

#############################
# 6. WEIGHTED REGRESSION
#############################

fit <- lm(y ~ x, data = pairs, weights = w)
summary(fit)

# Observed correlation
obs_r <- cor(pairs$x, pairs$y)
obs_r

#############################
# 6. GRID SEARCH
#############################
results <- expand.grid(
  radius = c(1, 2, 4, 8, 12),
  minx = c(1e-6, 1e-5, 1e-3),
  miny = c(1e-6, 1e-5, 1e-3),
  delfac = c(0.1, 0.5, 1, 2, 5)
)

grid_results <- apply(results, 1, function(par){
  
  radius <- par["radius"]
  minx   <- par["minx"]
  miny   <- par["miny"]
  delfac <- par["delfac"]
  
  # rodar find_pairs com os parâmetros do grid
  pairs <- find_pairs(tree, incX, incY,
                      radius = radius,
                      minx = minx,
                      miny = miny,
                      delfac = delfac)
  
  # remover NA
  pairs <- na.omit(pairs)
  
  # evitar casos sem dados suficientes
  if(nrow(pairs) < 5){
    return(data.frame(
      radius = radius,
      minx = minx,
      miny = miny,
      delfac = delfac,
      n_pairs = nrow(pairs),
      r = NA,
      slope = NA,
      p = NA
    ))
  }
  
  # modelo
  fit <- lm(y ~ x, data = pairs, weights = w)
  
  # correlação
  r <- cor(pairs$x, pairs$y)
  
  # extrair métricas
  coef_fit <- summary(fit)$coefficients
  
  data.frame(
    radius = radius,
    minx = minx,
    miny = miny,
    delfac = delfac,
    n_pairs = nrow(pairs),
    r = r,
    slope = coef_fit["x", "Estimate"],
    p = coef_fit["x", "Pr(>|t|)"],
    r2 = summary(fit)$r.squared
  )
})

grid_results <- bind_rows(grid_results)

head(grid_results)
summary(grid_results$r)
summary(grid_results$slope)

#############################
# 7. PERMUTATION TEST
#############################

random_test <- function(tree, incx, incy, reps = 1000){
  
  pairs_obs <- find_pairs(tree, incx, incy)
  pairs_obs <- na.omit(pairs_obs)
  
  obs_r <- cor(pairs_obs$x, pairs_obs$y)
  
  rand_r <- numeric(reps)
  
  for(i in 1:reps){
    
    incy_rand <- sample(incy)
    
    pairs_rand <- find_pairs(tree, incx, incy_rand)
    pairs_rand <- na.omit(pairs_rand)
    
    rand_r[i] <- cor(pairs_rand$x, pairs_rand$y)
    
  }
  
  p_value <- mean(abs(rand_r) >= abs(obs_r))
  
  return(list(
    obs_r = obs_r,
    p = p_value,
    rand_distribution = rand_r
  ))
}

#############################
# 8. RUN TEST
#############################

res <- random_test(tree, incX, incY, reps = 1000)

#############################
# 9. PLOT
#############################

plot(pairs$x, pairs$y,
     xlab = "Change in X",
     ylab = "Change in Y",
     main = "Lagged evolutionary correlation")

abline(fit, col = "black", lwd = 2)
abline(h = 0, v = 0, lty = 2, col = "gray")

#############################
# 10. NULL DISTRIBUTION
#############################

hist(res$rand_distribution,
     main = "Null distribution of correlations",
     xlab = "Correlation",
     xlim = range(c(res$rand_distribution, res$obs_r)))

abline(v = res$obs_r, col = "red", lwd = 3)

############################################################
# END
############################################################
# PLOT HAPLORRHINI
fit <- lm(y ~ x, data = pairs)
r2 <- summary(fit)$r.squared

p <- ggplot(pairs, aes(x = x, y = y)) +
  geom_point(size = 4, alpha = 0.8, color = "black") +
  
  # linha de regressão
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  
  # linhas de referência
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  labs(
    x = "Change in NMDS2",
    y = "Change in Sexual Dimorphism") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14)
  )

p = p +
  annotate("text",
           x = Inf, y = -Inf,
           label = paste0("R² = ", round(r2, 3)),
           hjust = 1.1, vjust = -0.5,
           size = 5)

p
# Salva o gráfico em alta resolução
ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/DELCOR_ScatterPlot_nmds2.png",
  plot = p,
  width = 12,
  height = 7,
  dpi = 300
)

df <- data.frame(rand = res$rand_distribution)
p <- ggplot(df, aes(x = rand)) +
  geom_histogram(bins = 30, fill = "gray", color = "white") +
  
  # linha do valor observado
  geom_vline(xintercept = res$obs_r, color = "red", linewidth = 1.5) +
  
  labs(
    x = "Correlation",
    y = "Frequency"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14)
  )

p = p + annotate("text",
             x = res$obs_r + 0.08,  # desloca 0.02 unidades para a direita
             y = Inf,
             label = "Observed",
             vjust = 2,
             color = "red",
             size = 6)

p
# Salva o gráfico em alta resolução
ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/matings/DELCOR_Histogram_nmds2.png",
  plot = p,
  width = 12,
  height = 7,
  dpi = 300
)
