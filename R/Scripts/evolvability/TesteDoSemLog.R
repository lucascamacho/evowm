# Different error magnitudes in evolvability measures
# final plots: 6 plots of evolvability x sexual dimorphism with different error sizes
#
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

library(ggplot2)
library(ggpmisc)
library(patchwork)

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}
prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# Read VCV matrices
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

# get species
matings = readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/Haplorrhini_MDS_Matings.RDS")

# Define species based on match between vcv and matings
commom <- intersect(names(vcv), matings$genus)


# Different errors
sigmas <- seq(0.01, 0.4, length.out = 9)

medidas <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/matings/averages_PCS_autovalues_primates.RDS")
species <- medidas$Genus

#
sp_comuns <- Reduce(intersect, list(
  names(vcv),
  medidas$Genus
))

results <- list()  # lista vazia
counter <- 1
for(i in 1:length(sigmas)){
  error <- sigmas[i]
  
  for(j in 1:length(sp_comuns)){    
    # choose species and vcv
    sp <- sp_comuns[j]
    
    # vcv
    covar <- vcv[[sp]]
    covar <- as.matrix(covar)
    
    # zero first eingenvector?
    #eig <- eigen(covar)
    #D <- diag(eig$values)
    #V <- eig$vectors
    #D2 <- D
    #D2[1,1]<- 0
    #covar <- V %*% D2 %*% t(V)
    
    # or remove it?
    #eig <- eigen(covar)
    #values <- eig$values[-1]
    #vectors <- eig$vectors[, -1]
    #covar <- vectors %*% diag(values) %*% t(vectors)
    
    # medias
    medias <- medidas$ByTrait_Averages[[sp]]
    
    #sd
    size <- (geomean(medias$Machos) + geomean( medias$Fêmeas)) / 2
    sd <- (medias$Machos - medias$Fêmeas) / size
    
    # ratio
    rat = geomean(medias$Fêmeas) / geomean(medias$Machos)
    
    # size matters the most
    x <- (medias$Machos[names(medias$Fêmeas)] + medias$Fêmeas) / 2
    D <- diag(x)
    P_scaled <- solve(D) %*% covar %*% solve(D)
    covar <- P_scaled
    
    # evolvability
    dimor_norm <- sd / sqrt(sum(sd ^ 2))
    evolv_normal <- as.numeric(t(dimor_norm) %*% covar %*% dimor_norm)
    
    # evolvability with error sigma[i]
    P <- covar
    B <- dimor_norm
    e <- rnorm(length(sd), sd = error)
    
    B_err <- B + e
    
    evolv_error <- as.numeric(
      t(B_err) %*% P %*% B_err /
        (t(B_err) %*% B_err)
    )
    
    #
    results[[counter]] <- data.frame(
      Species = sp,
      Sexual_Dimorphism = abs(norma(sd)) * sign(rat - 1),
      Evolvability = evolv_normal,
      Error_Evolvability = evolv_error,
      Sigma_Value = error,
      check.names = FALSE
    )
    counter <- counter + 1
  }
}

evolvas_error <- do.call(rbind, results)
rownames(evolvas_error) <- NULL

colnames(evolvas_error) <- c("Species", "Sexual_Dimorphism", "Evolvability", "Error_Evolvability",
                             "Sigma_Value")

#saveRDS(evolvas_error, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Errors.RDS")
#evolvas_error <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Errors.RDS")

# plot
p1 <- ggplot(evolvas_error, aes(x = Sexual_Dimorphism, y = Error_Evolvability)) +
   geom_point() +
   facet_wrap(~ Sigma_Value) +
   stat_poly_eq(
     aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
     formula = y ~ x,
     parse = TRUE,
     size = 5,
     label.x = "right"
   ) +
   labs(
     x = "Sexual Dimorphism",
     y = "Evolvability Estimation",
   ) +
   geom_smooth(method = "lm", se = TRUE) +
   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
   theme_classic(base_size = 14) +
   theme(
     panel.border = element_rect(color = "black", fill = NA),
     legend.position = "none",
     panel.grid = element_blank(),
     plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
     axis.title = element_text(size = 16, face = "bold"),
     axis.title.x = element_text(size = 16, face = "bold"),
     axis.title.y = element_text(size = 16, face = "bold"),
     axis.text = element_text(size = 14)
   )

p1

# ggsave(
#   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Sygmas_Different_Values.png",
#   plot = p1,
#   width = 12,
#   height = 7,
#   dpi = 300
# )

# # right names
evolvas_error$Species <- gsub("_", " ", evolvas_error$Species)

p2 <- ggplot(evolvas_error, aes(x = Species, y = Evolvability)) +
   geom_col() +
   scale_y_continuous(expand = c(0, 0)) +
   scale_x_discrete(drop = TRUE) +
   xlab(" Simian Species") +
   ylab("Evolvability of Sexual Dimorphism") +
   coord_flip() +
   theme_classic(base_size = 14) +
   theme(
     panel.border = element_rect(color = "black", fill = NA),
     legend.position = "none",
     panel.grid = element_blank(),
     plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
     axis.title = element_text(size = 16, face = "bold"),
     axis.title.x = element_text(size = 16, face = "bold"),
     axis.title.y = element_text(size = 16, face = "bold"),
     axis.text = element_text(size = 10)
   )

p2

# ggsave(
#   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Species.png",
#   plot = p2,
#   width = 12,
#   height = 7,
#   dpi = 300
# )

p3 <- ggplot(data = evolvas_error) +
  #geom_point(aes(x = Sexual_Dimorphism, y = Evolvability), size = 2) +
  geom_text(aes(x = Sexual_Dimorphism, y = Evolvability, label = Species), size = 2) +
  stat_poly_eq(
    aes(x = Sexual_Dimorphism, y = Evolvability, label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 5,
    label.x = "right"
  ) +
  geom_smooth(aes(x = Sexual_Dimorphism, y = Evolvability), method = "lm", se = TRUE) +
  xlab("Sexual Dimorfism") +
  ylab("Evolvability Estimates of Sexual Dimorphism") +
  theme_classic() +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )

p_combined <- p3 + p2
p_combined

# ggsave(
#    "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/SemLog.pdf",
#    plot = p_combined,
#    width = 18,
#    height = 10,
#    dpi = 300
#  )
