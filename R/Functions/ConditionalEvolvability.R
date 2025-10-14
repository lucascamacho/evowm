ConditionalEvolvability <- function(cov.matrix, beta.mat = NULL, iterations = 1000, size_skull) {
  # Funções auxiliares internas
  Norm <- function(x) {
    sqrt(sum(x * x))
  }
  
  Normalize <- function(x) {
    x / Norm(x)
  }
  
  # Converte para matriz da classe Matrix (para nearPD)
  cov.matrix <- Matrix::Matrix(cov.matrix)
  num.traits <- dim(cov.matrix)[1]
  
  # Se beta.mat não for fornecida, sorteia vetores aleatórios e normaliza
  if (is.null(beta.mat)) {
    beta.mat <- array(rnorm(num.traits * iterations), c(num.traits, iterations))
    beta.mat <- apply(beta.mat, 2, Normalize)
  }
  
  # Checa se a matriz é positiva definida, ajusta se necessário
  cov.matrix <- tryCatch({
    chol(cov.matrix)
    cov.matrix
  }, error = function(cond) {
    warning("Matrix is singular, using nearPD; results may be approximate")
    cov.matrix <- Matrix::nearPD(cov.matrix)[[1]]
    chol(cov.matrix)
  })
  
  # Calcula evolvabilidade condicional
  respostas <- 1 / diag(t(beta.mat) %*% solve(cov.matrix, beta.mat))
  
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