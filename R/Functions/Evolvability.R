Evolvability <- function(cov.matrix, beta.mat = NULL, iterations = 1000, size_skull) {
  # Funções auxiliares internas
  Norm <- function(x) {
    sqrt(sum(x * x))
  }
  
  Normalize <- function(x) {
    x / Norm(x)
  }
  
  # Número de traços
  num.traits <- dim(cov.matrix)[1]
  
  # Se beta.mat não for fornecida, sorteia vetores aleatórios e normaliza
  if (is.null(beta.mat)) {
    beta.mat <- array(rnorm(num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply(beta.mat, 2, Normalize)
  }
  
  # Calcula evolvabilidade (resposta)
  respostas <- diag(t(beta.mat) %*% cov.matrix %*% beta.mat)
  
  # Métricas
  respostas_pure <- mean(respostas)
  respostas_normal <- mean(respostas / size_skull)
  icv <- sd(respostas) / mean(respostas)
  icv_normal <- icv / size_skull
  
  # Retorna lista com resultados
  return(list(
    respostas = respostas,
    respostas_normal = respostas_normal,
    icv = icv,
    icv_normal = icv_normal,
    respostas_pure = respostas_pure
  ))
}