get_var_explained <- function(covar_list) {
  # Para cada matriz da lista
  res <- lapply(names(covar_list), function(species_name) {
    
    M <- covar_list[[species_name]]
    eig <- eigen(M, symmetric = TRUE)
    
    # Ordenar autovalores do maior para o menor
    eigvals <- eig$values
    eigvals_sorted <- sort(eigvals, decreasing = TRUE)
    
    # Proporção de variância explicada
    var_exp <- eigvals_sorted / sum(eigvals_sorted)
    
    # Criar dataframe para a espécie
    data.frame(
      species = species_name,
      rank = seq_along(var_exp),
      var_explained = var_exp,
      eigenvalue = eigvals_sorted
    )
  })
  
  # Empilha tudo
  do.call(rbind, res)
}
# temporario
df_var_exp = get_var_explained(vcv)

summary(df_var_exp$var_explained)

p1 = ggplot(df_var_exp, aes(x = rank, y = var_explained, color = species)) +
  geom_line(alpha = 0.7) +
  labs(y = "Proportion of variance explained", x = "Eigenvector rank") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14)
  )

p1

#ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Prop_Var_Explained.png", plot = p1,
#       width = 12,    # largura em inches
#       height = 8,   # altura em inches
#       dpi = 200)    # resolução

p1 = ggplot(df_var_exp, aes(x = rank, y = log(var_explained))) +
  geom_point(alpha = 0.1) +
  labs(y = "Natural Log of the Proportion of variance explained", x = "Eigenvector rank") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.4, alpha = 0.2) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14)
  )

p1

#summary(lm(log(df_var_exp$var_explained) ~ df_var_exp$rank))

#ggsave("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Log_Prop_Var_Explained.png", plot = p1,
#       width = 12,    # largura em inches
#       height = 8,   # altura em inches
#       dpi = 200)    # resolução
