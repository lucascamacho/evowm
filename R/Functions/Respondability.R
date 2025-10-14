# Calculate Respondability from a covariance matrix 
# size_skull is the average size of the populations skull
#
Respondability <- function(cov.matrix, beta.mat = NULL, iterations = 1000, size_skull) {
  # Intern functions
  Norm <- function(x) {
    sqrt(sum(x * x))
  }
  
  Normalize <- function(x) {
    x / Norm(x)
  }
  
  # N traits
  num.traits <- dim(cov.matrix)[1]
  
  # Se beta.mat não for fornecida, sorteia vetores aleatórios e normaliza
  if (is.null(beta.mat)) {
    beta.mat <- array(rnorm(num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply(beta.mat, 2, Normalize)
  }
  
  # Calcula resposta evolutiva
  respostas <- apply(cov.matrix %*% beta.mat, 2, Norm)
  
  # Métricas de respondabilidade
  respostas_pure <- mean(respostas)
  respostas_normal <- mean(respostas / size_skull)
  icv <- sd(respostas) / mean(respostas)
  icv_normal <- icv / size_skull
  
  # Retorna lista com os resultados
  return(list(
    respostas = respostas,
    respostas_pure = respostas_pure,
    respostas_normal = respostas_normal,
    icv = icv,
    icv_normal = icv_normal
  ))
}