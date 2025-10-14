# Função para checar uma única matriz P
check_matrix <- function(P) {
  P <- as.matrix(P)
  ev <- eigen(P, symmetric = TRUE)$values
  
  symmetric <- isTRUE(all.equal(P, t(P), tolerance = 1e-8))
  positive_definite <- all(ev > 0)
  chol_ok <- tryCatch({ chol(P); TRUE }, error = function(e) FALSE)
  
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