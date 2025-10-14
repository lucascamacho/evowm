Bias <- function(cov.matrix, beta.mat = NULL, iterations = 1000) {
  # Funções auxiliares internas
  Norm <- function(x) {
    sqrt(sum(x * x))
  }
  
  Normalize <- function(x) {
    x / Norm(x)
  }
  
  prod_interno <- function(x, y) {
    sum(x * y)
  }
  
  norma <- function(x) {
    sqrt(prod_interno(x, x))
  }
  
  corVector <- function(x, y) {
    prod_interno(x, y) / (norma(x) * norma(y))
  }
  
  # Número de traços
  num.traits <- dim(cov.matrix)[1]
  
  # Se beta.mat não for fornecida, sorteia vetores aleatórios e normaliza
  if (is.null(beta.mat)) {
    beta.mat <- array(rnorm(num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply(beta.mat, 2, Normalize)
  }
  
  # Vetor para armazenar correlações
  cor_v <- numeric(iterations)
  
  # Calcula correlação de resposta com o primeiro autovetor da matriz
  for (k in 1:iterations) {
    resposta <- cov.matrix %*% beta.mat[, k]
    resposta <- Normalize(resposta)
    cor_v[k] <- abs(corVector(resposta, eigen(cov.matrix)$vectors[, 1]))
  }
  
  # Métricas
  med_respo <- mean(cor_v)
  icv <- sd(cor_v) / mean(cor_v)
  
  # Retorna lista com resultados
  return(list(
    cor_v = cor_v,
    med_respo = med_respo,
    icv = icv
  ))
}