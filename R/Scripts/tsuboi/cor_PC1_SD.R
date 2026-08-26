# Align SD PCmax
# vector correlation between sexual dimorphism vector and the PC1 of the average P
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

library(ape)
library(evolqg)
library(ggplot2)
library(ggrepel)

load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/vcv.RData")

load("~/Desktop/Primatrees.RData")
species <- names(vcv)
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

medias <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData")
sd <- load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData")

clades <- read.csv("~/Desktop/taxonomy.csv")

geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

all_cov_matrices <- PhyloW(tree, vcv)
ancestral <- getMRCA(tree, tree$tip.label)
anc_vcv <- all_cov_matrices[[as.character(ancestral)]]

eig <- eigen(anc_vcv)
D <- diag(eig$values)
V <- eig$vectors
D2 <- D
D2[1,1] <- max(eig$values) * 1e-8
anc_vcv <- V %*% D2 %*% t(V)

pcmax <- eigen(anc_vcv)$vectors[,1]

# Extract the 39 sexual dimorphism variables
sd_vectors <- sexual_dimorphism_df[, as.character(1:39)]

# Calculate vector correlation and norm for each genus
results_sd_pcmax <- data.frame(
  Genus = sexual_dimorphism_df$Genus,
  Vector_Correlation = apply(
    sd_vectors,
    1,
    function(x) corVector(x, pcmax)
  ),
  Norm = apply(
    sd_vectors,
    1,
    norma
  )
)

# Add clade information
results_sd_pcmax <- results_sd_pcmax %>%
  left_join(
    clades %>%
      select(GENUS, CLADE) %>%
      distinct(),
    by = c("Genus" = "GENUS")
  )

p_sd_pcmax <- ggplot(
  results_sd_pcmax,
  aes(
    x = atanh(Vector_Correlation),
    y = Norm,
    color = CLADE
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 3
  ) +
  
  geom_text_repel(
    aes(label = Genus),
    size = 3,
    show.legend = FALSE
  ) +
  
  labs(
    x = "Z transformed vector correlation of log SD with average PCmax (size removed)",
    y = "Sexual dimorphism vector norm",
    color = "Clade"
  ) +
  
  theme_classic()

p_sd_pcmax

ggsave(
  "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/CorVector_PCmax_SD.png",
  plot = p_sd_pcmax,
  width = 12,
  height = 7,
  dpi = 300
)
