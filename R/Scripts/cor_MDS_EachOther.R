if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(reshape2)){install.packages("reshape2"); library(reshape2)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}
if(!require(purrr)){install.packages("purrr"); library(purrr)}

mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")

# CATARRHINI
# correlation between MDS and original variables (correlation matrix)
mds = mds[which(mds$dados.PARVORDER == "Catarrhini"),]

cor_mat <- sapply(mds[, 5:14], function(col) {
  c(
    cor(mds[[1]], col, use = "complete.obs"),
    cor(mds[[2]], col, use = "complete.obs")
  )
})

rownames(cor_mat) = c("NMDS1", "NMDS2")
colnames(cor_mat) = c("SOCIAL_ORG", "MATING_SYSTEM", "PROP_MALES_FEMALES", "DOMINANCE",
                      "AGGRESSION", "ALL_MALE_GROUPS", "FURTIVE_COPULATIONS", "INFANTICIDE",
                      "MULTILEVEL_SOCIETY", "NUMBER_PAPERS")


# p-value matrix (same structure, but with p-values)
p_mat <- sapply(mds[, 5:14], function(col) {
  c(
    cor.test(mds[[1]], col)$p.value,
    cor.test(mds[[2]], col)$p.value
  )
})

rownames(p_mat) = c("NMDS1", "NMDS2")
colnames(p_mat) = c("SOCIAL_ORG", "MATING_SYSTEM", "PROP_MALES_FEMALES", "DOMINANCE",
                    "AGGRESSION", "ALL_MALE_GROUPS", "FURTIVE_COPULATIONS", "INFANTICIDE",
                    "MULTILEVEL_SOCIETY", "NUMBER_PAPERS")

# Derrete as matrizes
cor_long <- melt(cor_mat)
p_long <- melt(p_mat)

# Junta as duas
df <- left_join(cor_long, p_long, by = c("Var1", "Var2"))
colnames(df) <- c("NMDS", "Variable", "Correlation", "p_value")

# Cria a coluna com o formato do texto
df <- df %>%
  mutate(label = sprintf("%.2f", Correlation),
         style = case_when(
           p_value < 0.01 ~ "bold",
           p_value < 0.05 ~ "italic",
           TRUE ~ "plain"
         ))

# --- Plot ---
c = ggplot(df, aes(x = Variable, y = NMDS, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label, fontface = style), size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Correlation") +
  ggtitle("Catarrhini") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title = element_blank(),
    legend.position = "none"
  )

# PLATYRRHINI
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
mds = mds[which(mds$dados.PARVORDER == "Platyrrhini"),]

cor_mat <- sapply(mds[, 5:14], function(col) {
  c(
    cor(mds[[1]], col, use = "complete.obs"),
    cor(mds[[2]], col, use = "complete.obs")
  )
})

rownames(cor_mat) = c("NMDS1", "NMDS2")
colnames(cor_mat) = c("SOCIAL_ORG", "MATING_SYSTEM", "PROP_MALES_FEMALES", "DOMINANCE",
                      "AGGRESSION", "ALL_MALE_GROUPS", "FURTIVE_COPULATIONS", "INFANTICIDE",
                      "MULTILEVEL_SOCIETY", "NUMBER_PAPERS")


# p-value matrix (same structure, but with p-values)
p_mat <- sapply(mds[, 5:14], function(col) {
  c(
    cor.test(mds[[1]], col)$p.value,
    cor.test(mds[[2]], col)$p.value
  )
})

rownames(p_mat) = c("NMDS1", "NMDS2")
colnames(p_mat) = c("SOCIAL_ORG", "MATING_SYSTEM", "PROP_MALES_FEMALES", "DOMINANCE",
                    "AGGRESSION", "ALL_MALE_GROUPS", "FURTIVE_COPULATIONS", "INFANTICIDE",
                    "MULTILEVEL_SOCIETY", "NUMBER_PAPERS")

# Derrete as matrizes
cor_long <- melt(cor_mat)
p_long <- melt(p_mat)

# Junta as duas
df <- left_join(cor_long, p_long, by = c("Var1", "Var2"))
colnames(df) <- c("NMDS", "Variable", "Correlation", "p_value")

# Cria a coluna com o formato do texto
df <- df %>%
  mutate(label = sprintf("%.2f", Correlation),
         style = case_when(
           p_value < 0.01 ~ "bold",
           p_value < 0.05 ~ "italic",
           TRUE ~ "plain"
         ))

# --- Plot ---
p = ggplot(df, aes(x = Variable, y = NMDS, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label, fontface = style), size = 5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Correlation") +
  ggtitle("Platyrrhini") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.title = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

fig_final = c + p
fig_final

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_MDS_Original_Variables.png", plot = fig_final,
       width = 18,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução


# correlation between original variables
# CATARRHINI
subdata <- mds[1:38, 5:14]

# calcular a matriz de correlação
cor_mat <- cor(subdata, use = "complete.obs")
colnames(cor_mat) <- colnames(cor_mat) %>%
  gsub("^dados\\.", "", .) %>%  # remove "dados." do início
  gsub("_", " ", .)             # troca "_" por espaço

rownames(cor_mat) <- rownames(cor_mat) %>%
  gsub("^dados\\.", "", .) %>%  # remove "dados." do início
  gsub("_", " ", .)             # troca "_" por espaço

vars_data <- mds[1:38, 5:14]   # ajustar colunas conforme necessário

# nomes limpos
vars <- colnames(vars_data) %>%
  gsub("^dados\\.", "", .) %>%
  gsub("_", " ", .)

n <- ncol(vars_data)

# Matriz de correlação
cor_mat <- cor(vars_data, use = "pairwise.complete.obs")
colnames(cor_mat) <- vars
rownames(cor_mat) <- vars

# Matriz de p-values
p_mat <- matrix(NA, nrow = n, ncol = n)
colnames(p_mat) <- vars
rownames(p_mat) <- vars

for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      p_mat[i, j] <- NA
    } else {
      test <- cor.test(vars_data[, i], vars_data[, j])
      p_mat[i, j] <- test$p.value
    }
  }
}

# Transformar em long format
cor_df <- melt(cor_mat, varnames = c("var1", "var2"), value.name = "correlacao")
p_df   <- melt(p_mat, varnames = c("var1", "var2"), value.name = "p_value")


# criar plot_df juntando correlações e p-values
plot_df <- cor_df %>%
  left_join(p_df, by = c("var1", "var2")) %>%
  # filtrar para upper triangle se quiser só metade
  filter(match(var1, vars) <= match(var2, vars)) %>%
  mutate(
    label_text = round(correlacao, 2),
    font_style = case_when(
      !is.na(p_value) & p_value < 0.01 ~ "bold",
      !is.na(p_value) & p_value >= 0.01 & p_value < 0.05 ~ "italic",
      TRUE ~ "plain"
    )
  )

# Plot
c <- ggplot(plot_df, aes(x = var1, y = var2, fill = correlacao)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label_text, fontface = font_style), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "white") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Catarrhini",
    x = "",
    y = "",
    fill = "Correlation"
  ) +
  scale_y_discrete(limits = rev(vars))   # inverter eixo Y

# Platyrrhini
subdata <- mds[39:62, 5:14]

# calcular a matriz de correlação
cor_mat <- cor(subdata, use = "complete.obs")
colnames(cor_mat) <- colnames(cor_mat) %>%
  gsub("^dados\\.", "", .) %>%  # remove "dados." do início
  gsub("_", " ", .)             # troca "_" por espaço

rownames(cor_mat) <- rownames(cor_mat) %>%
  gsub("^dados\\.", "", .) %>%  # remove "dados." do início
  gsub("_", " ", .)             # troca "_" por espaço

vars_data <- mds[39:62, 5:14]   # ajustar colunas conforme necessário

# nomes limpos
vars <- colnames(vars_data) %>%
  gsub("^dados\\.", "", .) %>%
  gsub("_", " ", .)

n <- ncol(vars_data)

# Matriz de correlação
cor_mat <- cor(vars_data, use = "pairwise.complete.obs")
colnames(cor_mat) <- vars
rownames(cor_mat) <- vars

# Matriz de p-values
p_mat <- matrix(NA, nrow = n, ncol = n)
colnames(p_mat) <- vars
rownames(p_mat) <- vars

for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      p_mat[i, j] <- NA
    } else {
      test <- cor.test(vars_data[, i], vars_data[, j])
      p_mat[i, j] <- test$p.value
    }
  }
}

# Transformar em long format
cor_df <- melt(cor_mat, varnames = c("var1", "var2"), value.name = "correlacao")
p_df   <- melt(p_mat, varnames = c("var1", "var2"), value.name = "p_value")


# criar plot_df juntando correlações e p-values
plot_df <- cor_df %>%
  left_join(p_df, by = c("var1", "var2")) %>%
  # filtrar para upper triangle se quiser só metade
  filter(match(var1, vars) <= match(var2, vars)) %>%
  mutate(
    label_text = round(correlacao, 2),
    font_style = case_when(
      !is.na(p_value) & p_value < 0.01 ~ "bold",
      !is.na(p_value) & p_value >= 0.01 & p_value < 0.05 ~ "italic",
      TRUE ~ "plain"
    )
  )

# Plot
p <- ggplot(plot_df, aes(x = var1, y = var2, fill = correlacao)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label_text, fontface = font_style), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "white") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        #axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Platyrrhini",
    x = "",
    y = "",
    fill = "Correlation"
  ) +
  scale_y_discrete(limits = rev(vars))   # inverter eixo Y


fig_final = c + p
fig_final

ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_MDS_EachOther.png", plot = fig_final,
       width = 18,    # largura em inches
       height = 8,   # altura em inches
       dpi = 200)    # resolução
