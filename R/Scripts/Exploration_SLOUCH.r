# CATARRHINI
load("~/Dropbox/Doc/Code/evowm/R/Outputs/catarrhini_fits.RData")

library(ggplot2)
library(ape)
library(dplyr)
library(tidyr)
library(phytools)
library(scales)

# INTEGRATION
evo_slope <- 0.5009472
se <- 0.1244806
intercept <- 0.03858323  # vindo do summary do modelo

lm_fit <- lm(response ~ integration, data = df_ou)

# Extrai coeficientes e R²
intercept_val <- coef(lm_fit)[1]
slope_val <- coef(lm_fit)[2]
r2_val <- summary(lm_fit)$r.squared

d <- ggplot(df_ou, aes(x = integration, y = response)) +
  geom_point(size = 2.5, color = "#1b9e77", alpha = 0.9, stroke = 2) +
  geom_abline(
    intercept = intercept_val,
    slope = slope_val,
    size = 1.2,
    color = "#1b9e77"
  ) +
  geom_abline(
    intercept = intercept_val,
    slope = slope_val + 1.96 * se,
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_abline(
    intercept = intercept_val,
    slope = slope_val - 1.96 * se,
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  
  # Adiciona texto com os coeficientes e R²
  annotate(
    "text",
    x = Inf, y = -Inf,                      # canto inferior direito
    hjust = 1.1, vjust = -0.5,
    label = paste0(
      "Slope = ", round(slope_val, 3), "\n",
      "Intercept ≈ ", round(intercept_val, 3), "\n",
      " R² = ", round(r2_val, 3)
    ),
    parse = FALSE,
    size = 5,
    color = "black"
  ) +
  
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14)
  ) +
  labs(
    x = "Integration index",
    y = "Average Normalized Sexual Dimorphism"
  )

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Integration_SD_OU.png", plot = d,
       width = 16,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução

### OU CURVE
alpha <- 0.05
sigma2 <- 0.05072902
stationary_var <- 0.5072902
half_life <- 13.86294
tree_height <- max(branching.times(tree_c))

# Tempo para plotar (por exemplo, até 50 unidades)
t <- rev(seq(0, tree_height, by = 0.1))

# Variância ao longo do tempo
var_t <- (sigma2 / (2 * alpha)) * (1 - exp(-2 * alpha * t))

# Data frame para plotar
df_var <- data.frame(
  time = t,
  variance = var_t
)

prop_half_life <- half_life / tree_height
prop_half_life

h <- ggplot(df_var, aes(x = time, y = variance)) +
  geom_line(color = "#1b9e77", size = 1.2) +
  geom_hline(yintercept = stationary_var, linetype = "dashed", color = "red") +
  geom_vline(xintercept = half_life, linetype = "dotted", color = "blue") +
  annotate(
    "text",
    x = half_life + 2,
    y = stationary_var * 0.5,
    label = paste0("Half-life ≈ ", round(half_life, 2)),
    color = "blue",
    hjust = 0,
    size = 6   # 🔹 aumenta o tamanho da fonte
  ) +
  annotate(
    "text",
    x = max(t) * 0.7,
    y = stationary_var + 0.02,
    label = paste0("Stationary Variance ≈ ", round(stationary_var, 3)),
    color = "red",
    size = 6   # 🔹 aumenta o tamanho da fonte
  ) +
  labs(
    x = "Evolutionary Time in Millions of Years",
    y = "Cumulative Variance",
    title = "Variance Accumulation Through Evolutionary Time"  # 🔹 adiciona título
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5) # 🔹 centraliza e aumenta o título
  )  

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolution_Curve_OU.png", plot = h,
       width = 16,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


####
# --- Criar os data frames com os slopes ---
opt <- data.frame(
  predictor = "integration",
  estimate  = 1.047039,
  se        = 0.260179,
  type      = "Optimal"
)

evo <- data.frame(
  predictor = "integration",
  estimate  = 0.5009472,
  se        = 0.1244806,
  type      = "Evolutionary"
)

# juntar
df <- bind_rows(opt, evo)

# calcular IC 95% aproximado e se IC inclui 0
df <- df %>%
  mutate(
    ci_low  = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se,
    sig     = ifelse(ci_low > 0 | ci_high < 0, TRUE, FALSE)
  )

# Ajustar posição do ★ de forma proporcional
df <- df %>%
  mutate(star_y = ci_high + 0.05 * (max(ci_high) - min(ci_low)))

p = ggplot(df, aes(x = type, y = estimate, color = type)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.1, size = 1) +
  geom_text(aes(y = star_y, label = ifelse(sig, "★", "")), 
            vjust = -0.5, size = 6, show.legend = FALSE) +
  labs(
    x = "",
    y = "Slope value",
    title = "Optimal vs Evolutionary slope for Integration",
    subtitle = "Error bars = 95% CI (±1.96·SE). ★ = CI does not include 0"
  ) +
  scale_color_manual(values = c("#1b9e77", "#1b9e77")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 16)
  )

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Optimum_Evolutive_Slope_OU.png", plot = p,
       width = 16,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


####
phy <- fits$OU_INTEGRATION$tree$phy
resid <- fits$OU_INTEGRATION$beta_evolutionary$residuals[,1]
names(resid) <- phy$tip.label

# Escala os resíduos para 1-100
res_scaled <- rescale(resid, to=c(1,100))

# Cria paleta azul-branco-vermelho
cols <- colorRampPalette(c("blue","white","red"))(100)

# Plota a árvore
plotTree(phy, fsize=0.8)
tiplabels(pch=21, bg=cols[round(res_scaled)], cex=1.5)


# PLATYRRHINI
load("~/Dropbox/Doc/Code/evowm/R/Outputs/platyrrhini_fits.RData")
summary(fits$BM)

sigma2    <- fits$BM$evolpar$sigma2_y
intercept <- fits$BM$beta_primary$coefficients[1]
residuals <- fits$BM$beta_primary$residuals[,1]
tree      <- fits$BM$tree$phy

species   <- fits$BM$tree$phy$tip.label

# -------------------------------
# 2️⃣ Plotar resíduos na árvore
# -------------------------------
# Cria um vetor nomeado com os resíduos
res_vec <- setNames(residuals, species)

# Plot usando contMap para visualização
cont_map <- contMap(tree, res_vec, plot = FALSE)
plot(cont_map, fsize = 0.7, lwd = 3, legend = TRUE, outline = FALSE)

# -------------------------------
# 3️⃣ Histogramas de resíduos
# -------------------------------
res_df <- data.frame(residuals = residuals)

ggplot(res_df, aes(x = residuals)) +
  geom_histogram(bins = 15, fill = "#7570b3", color = "black", alpha = 0.8) +
  labs(
    x = "Residuals (BM)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(size = 18),     # nome do eixo X
        axis.title.y = element_text(size = 18))    # nome do eixo Y))

# -------------------------------
# 4️⃣ Simulações sob BM
# -------------------------------
# Simula 100 trajetórias sob BM
nsim <- 100
sim_data <- replicate(nsim, fastBM(tree, sig2 = sigma2, a = intercept))

# Converte para data.frame para ggplot
sim_df <- as.data.frame(sim_data)
sim_df$species <- rownames(sim_df)
sim_long <- sim_df %>% tidyr::pivot_longer(cols = -species, names_to = "sim", values_to = "trait")

ç = ggplot(sim_long, aes(x = trait)) +
  geom_histogram(bins = 20, fill = "#7570b3", alpha = 0.5, color = "black") +
  labs(
    x = "Simulated values (BM)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(size = 18),     # nome do eixo X
        axis.title.y = element_text(size = 18))    # nome do eixo Y))

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Simulated_BM_Values.png", plot = ç,
       width = 20,    # largura em inches
       height = 12,   # altura em inches
       dpi = 200)    # resolução

# -------------------------------
# 5️⃣ Observado vs. esperado
# -------------------------------
obs_vals <- intercept + residuals
obs_df <- data.frame(
  species = species,
  observed = obs_vals
)

m = ggplot(obs_df, aes(x = species, y = observed)) +
      geom_col(fill = "#7570b3", color = "black", alpha = 0.8) +
        labs(
        x = "Species",
        y = "Observed values (Intercept + residuals)") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 18),     # nome do eixo X
        axis.title.y = element_text(size = 18))    # nome do eixo Y))

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Observed_BM_Model.png", plot = m,
       width = 20,    # largura em inches
       height = 12,   # altura em inches
       dpi = 200)    # resolução

