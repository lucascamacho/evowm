# Different error magnitudes in evolvability measures
# final plots: 6 plots of evolvability x sexual dimorphism with different error sizes
#
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability")

library(ggplot2)
library(ggpmisc)
library(patchwork)

prod_interno = function(x, y) sum(x * y)
norma = function(x) sqrt(prod_interno(x, x))
corVector = function(x, y) prod_interno(x, y)/(norma(x)*norma(y))

# geometric mean function
geomean = function(vector){
  g = exp(mean(log(vector)))
  return(g)
}

# Read VCV matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
temp <- list.files(pattern = "*.csv")
vcv <- lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
names(vcv) <- gsub(".csv", replacement= "", temp)

# Different errors
sigmas <- seq(0.1, 0.9, length.out = 3)

medidas <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/averages_PCS_autovalues_primates.RDS")
species <- medidas$Species

#
sp_comuns <- Reduce(intersect, list(
  names(vcv),
  medidas$Species
))

results <- list()  # lista vazia
counter <- 1
for(i in 1:length(sigmas)){
  error <- sigmas[i]
  
  for(j in 1:length(sp_comuns)){    
    # choose species and vcv
    sp <- sp_comuns[j]
    
    medias <- medidas$ByTrait_Averages[[sp]]
    
    
    # vcv
    covar <- vcv[[sp]]
    covar <- as.matrix(covar)
    
    # calcula médias dos traits
    medias_por_trait <- (medias$Machos + medias$Fêmeas) / 2
    medias_por_trait <- as.numeric(medias_por_trait)
    mean_prod <- outer(medias_por_trait, medias_por_trait, "*") 
    mat_scaled <- covar / mean_prod   # padroniza tudo
    diag(mat_scaled) <- diag(covar) / (medias_por_trait ^ 2)   
    covar <- mat_scaled
    
    
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
    sd <- medias$Machos - medias$Fêmeas
    
    # ratio
    rat = geomean(medias$Fêmeas) / geomean(medias$Machos)
    
    # evolvability
    dimor_norm <- sd / sqrt(sum(sd ^ 2))
    evolv_normal <- as.numeric(t(dimor_norm) %*% covar %*% dimor_norm)
    
    # evolvability with error sigma[i]
    # error can be a product of variances
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
      Sexual_Dimorphism = norma(sd),
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
  geom_point(size = 6, alpha = 0.8) +
  facet_wrap("Sigma_Value") +
  stat_poly_eq(
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 10,
    label.x = "right"
  ) +
  labs(
    title = "Evolvability and Sexual Dimorphism with Different Sygma Values",
    x = "Sexual Dimorphism",
    y = "Evolvability Estimation"
  ) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title = element_text(size = 28, face = "bold"),
    axis.text = element_text(size = 24, face = "bold", colour = "black"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24, face = "bold"),
    strip.text = element_text(size = 24, face = "bold")
  )

p1

 ggsave(
   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Sygmas_Different_Values.png",
   plot = p1,
   width = 14,
   height = 7,
   dpi = 300
 )

# right names
evolvas_error$Species <- gsub("_", " ", evolvas_error$Species)

#
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
  geom_vline(xintercept = 0, col = "green", linetype = "dashed") +
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
#   "~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability_Species_SD.png",
#   plot = p_combined,
#   width = 18,
#   height = 10,
#   dpi = 300
# )