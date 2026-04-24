library(ggplot2)

# Read VCV matrices
setwd("~/Dropbox/Doc/Code/evowm/R/Outputs/log/")
temp = list.files(pattern = "*.csv")
vcv = lapply(temp, read.csv, header = TRUE, dec = ".", sep = ' ', row.names = 1)
names(vcv)  = gsub(".csv", replacement= "", temp)

# Define functions
simulate_dimension <- function(p, sigmas, nrep = 100) {
  
  # B verdadeiro
  B <- rnorm(p)
  
  # P aleatória PD
  A <- matrix(rnorm(p^2), p, p)
  P <- crossprod(A)
  
  theta_real <- as.numeric(t(B) %*% P %*% B / (t(B) %*% B))
  
  res <- expand.grid(
    sigma = sigmas,
    rep = 1:nrep
  )
  
  res$theta_obs <- NA
  
  for (i in seq_len(nrow(res))) {
    
    e <- rnorm(length(B), sd = res$sigma[i])
    
    num_error <- as.numeric(
      t(B) %*% P %*% e +
        t(e) %*% P %*% B +
        t(e) %*% P %*% e
    )
    
    den_error <- as.numeric(
      t(B) %*% B +
        t(B) %*% e +
        t(e) %*% B +
        t(e) %*% e
    )
    
    res$theta_obs[i] <- theta_real * (num_error / den_error)
  }
  
  res$theta_real <- theta_real
  res$p <- p
  return(res)
}

make_modular_P <- function(p, n_modules, rho_within = 0.7, rho_between = 0.1) {
  
  module_size <- p / n_modules
  stopifnot(module_size == round(module_size))
  
  P <- matrix(rho_between, p, p)
  
  for (m in 1:n_modules) {
    idx <- ((m - 1) * module_size + 1):(m * module_size)
    P[idx, idx] <- rho_within
  }
  
  diag(P) <- 1
  
  # garantir positiva definida
  eigenvals <- eigen(P, symmetric = TRUE)$values
  if (min(eigenvals) <= 0) {
    P <- P + diag(abs(min(eigenvals)) + 0.01, p)
  }
  
  return(P)
}

simulate_modularity <- function(P, B, sigmas, nrep = 100) {
  
  theta_real <- as.numeric(t(B) %*% P %*% B / (t(B) %*% B))
  
  res <- expand.grid(
    sigma = sigmas,
    rep = 1:nrep
  )
  
  res$theta_obs <- NA
  
  for (i in seq_len(nrow(res))) {
    
    e <- rnorm(length(B), sd = res$sigma[i])
    
    num_error <- as.numeric(
      t(B) %*% P %*% e +
        t(e) %*% P %*% B +
        t(e) %*% P %*% e
    )
    
    den_error <- as.numeric(
      t(B) %*% B +
        t(B) %*% e +
        t(e) %*% B +
        t(e) %*% e
    )
    
    res$theta_obs[i] <- theta_real * (num_error / den_error)
  }
  
  res$theta_real <- theta_real
  return(res)
}

# Different errors
sigmas <- seq(0.01, 1, length.out = 30)

# Different scenarios of p (number of traits)
res_all <- rbind(
  simulate_dimension(p = 5,  sigmas),
  simulate_dimension(p = 20, sigmas),
  simulate_dimension(p = 50, sigmas)
)

ggplot(res_all, aes(x = sigma, y = theta_obs)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun.data = median_hilow, geom = "ribbon", alpha = 0.2) +
  facet_wrap(~ p, scales = "free_y") +
  geom_hline(aes(yintercept = theta_real),
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Magnitude of the error (σ)",
    y = "Estimative of evolvability",
    title = "Different sizes of P matrix (number of traits)"
  )

p <- 38
B <- rnorm(p)
sigmas <- seq(0.01, 1, length.out = 30)

P_random <- crossprod(matrix(rnorm(p^2), p))
#P_modular <- make_modular_P(p, n_modules = 4)
P_empirical <- as.matrix(vcv[[1]])
p <- ncol(P_empirical)

medidas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")
#B <- rnorm(p)
B <- log(medidas$ByTrait_Averages$Cercocebus_agilis$Machos) - log(medidas$ByTrait_Averages$Cercocebus_agilis$Fêmeas)
B = B[-20]

P_empirical <- (P_empirical + t(P_empirical)) / 2
eig <- eigen(P_empirical, symmetric = TRUE)
if (min(eig$values) <= 0) {
  P_empirical <- P_empirical + diag(abs(min(eig$values)) + 1e-6, p)
}

P_random <- crossprod(matrix(rnorm(p^2), p))

res_rand <- simulate_modularity(P_random, B, sigmas)
res_emp  <- simulate_modularity(P_empirical, B, sigmas)

res_rand$type <- "Random P"
res_emp$type  <- "Empirical P"

res <- rbind(res_rand, res_emp)

ggplot(res, aes(x = sigma, y = theta_obs, color = type)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun.data = median_hilow,
               geom = "ribbon", alpha = 0.15, aes(fill = type)) +
  facet_wrap(~ type, scales = "free_y") +
  geom_hline(aes(yintercept = theta_real),
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Magnitude of the error (σ)",
    y = "Estimation of evolvability",
    title = "Erro sob P aleatória vs modular"
  )
