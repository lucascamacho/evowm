setwd("~/Dropbox/Doc/Data/vcv")

library(Matrix)

# Função para checar uma matriz P
check_matrix <- function(P) {
  P = as.matrix(P)
  # Eigenvalues
  ev <- eigen(P, symmetric = TRUE)$values
  
  # Simetria
  symmetric <- isTRUE(all.equal(P, t(P), tolerance = 1e-8))
  
  # Positiva definida (todas eigen > 0)
  positive_definite <- all(ev > 0)
  
  # Cholesky check (alternativa)
  chol_ok <- tryCatch({ chol(P); TRUE }, error=function(e) FALSE)
  
  # Estatísticas
  min_ev <- min(ev)
  max_ev <- max(ev)
  n_pos <- sum(ev > 0)
  kappa <- if (all(ev > 0)) max_ev / min_ev else NA
  logdet <- if (all(ev > 0)) sum(log(ev)) else NA
  detP <- det(P)
  
  return(list(
    symmetric = symmetric,
    positive_definite = positive_definite,
    chol_ok = chol_ok,
    min_ev = min_ev,
    max_ev = max_ev,
    n_pos = n_pos,
    kappa = kappa,
    logdet = logdet,
    det = detP
  ))
}


# read P matrices
temp = list.files(pattern="*.csv")
vcv = lapply(temp, read.csv, header = FALSE, dec = ",")
names(vcv)  = gsub(".csv", replacement = "", temp)

# Loop sobre as 39 matrizes
results_list <- lapply(vcv, check_matrix)

# Converte lista em data.frame
results_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
  cbind(matrix_id = i, as.data.frame(results_list[[i]]))
}))

# Visualizar resultados
print(results_df)