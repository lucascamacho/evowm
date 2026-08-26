setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

# load packages
library(ape)
library(evolqg)

# load tree
load("~/Desktop/Primatrees.RData")
genus <- tree$tip.label
#filename <- "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
#tree <- read.nexus(filename)

# load data
load("~/Desktop/Primaset_RawData.RData")

species <- names(vcv)
tree <- drop.tip(tree, setdiff(tree$tip.label, species))

calculate_P_matrices <- function(info, traits) {
  
  genera <- unique(info$GEN)
  genera <- genera[!is.na(genera)]
  
  vcv <- list()
  
  for (genus in genera) {
    
    rows <- which(info$GEN == genus)
    
    info_genus <- info[rows, ]
    traits_genus <- traits[rows, , drop = FALSE]
    
    # Combine data
    data_genus <- data.frame(
      traits_genus,
      SPE = info_genus$SPE,
      SEX = info_genus$SEX
    )
    
    # Remove missing predictors
    complete <- complete.cases(data_genus[, c("SPE", "SEX")])
    data_genus <- data_genus[complete, ]
    
    # Convert predictors to factors
    data_genus$SPE <- factor(data_genus$SPE)
    data_genus$SEX <- factor(data_genus$SEX)
    
    # Basic requirements
    if (nrow(data_genus) <= 39) {
      warning(
        genus, ": too few complete specimens (N = ",
        nrow(data_genus), "). Skipping."
      )
      next
    }
    
    # Fit model
    model <- try(
      lm(
        as.matrix(log(data_genus[, 1:39])) ~ SPE * SEX,
        data = data_genus
      ),
      silent = TRUE
    )
    
    # Check whether model worked
    if (inherits(model, "try-error")) {
      warning(genus, ": model failed. Skipping.")
      next
    }
    
    # Calculate P matrix
    P <- try(
      CalculateMatrix(model),
      silent = TRUE
    )
    
    # Check whether CalculateMatrix worked
    if (inherits(P, "try-error")) {
      warning(genus, ": CalculateMatrix failed. Skipping.")
      next
    }
    
    # Store matrix
    vcv[[genus]] <- P
    
    cat(
      "Finished:", genus,
      "| N =", nrow(data_genus),
      "| species =", nlevels(data_genus$SPE),
      "\n"
    )
  }
  
  return(vcv)
}

vcv <- calculate_P_matrices(info, traits)

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

#
check <- lapply(vcv, check_matrix)
do.call(rbind.data.frame, check)

problematic <- c(
  "Avahi",
  "Callicebus",
  "Brachyteles",
  "Cheracebus",
  "Allochrocebus"
)

for (g in problematic) {
  
  cat("\n\n==========", g, "==========\n")
  
  x <- info[info$GEN == g, ]
  
  print(table(x$SPE, x$SEX))
  
}

summary_models <- data.frame(
  GEN = character(),
  N = numeric(),
  N_species = numeric(),
  df_residual = numeric()
)

for (g in names(vcv)) {
  
  rows <- which(info$GEN == g)
  
  x <- data.frame(
    traits[rows, ],
    SPE = info$SPE[rows],
    SEX = info$SEX[rows]
  )
  
  x <- x[complete.cases(x[, c("SPE", "SEX")]), ]
  
  x$SPE <- factor(x$SPE)
  x$SEX <- factor(x$SEX)
  
  model <- lm(
    as.matrix(x[, 1:39]) ~ SPE * SEX,
    data = x
  )
  
  summary_models <- rbind(
    summary_models,
    data.frame(
      GEN = g,
      N = nrow(x),
      N_species = nlevels(x$SPE),
      df_residual = df.residual(model)
    )
  )
}

summary_models

# Removendo os generos problematic <- c(
#"Avahi",
#"Callicebus",
#"Brachyteles",
#"Cheracebus",
#"Allochrocebus"
#

problematic <- c(
  "Avahi",
  "Callicebus",
  "Brachyteles",
  "Cheracebus",
  "Allochrocebus"
)

vcv <- vcv[!names(vcv) %in% problematic]
save(vcv, file = "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/vcv.RData")
