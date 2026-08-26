
load("~/Desktop/Primaset_RawData.RData")

clades <- read.csv("~/Desktop/taxonomy.csv")

measure_info <- data.frame(
  Measure = distn,
  Module = c(
    "Face",       # IS-PM
    "Face",       # IS-NSL
    "Face",       # IS-PNS
    "Face",       # PM-ZS
    "Face",       # PM-ZI
    "Face",       # PM-MT
    "Face",       # NSL-NA
    "Face",       # NSL-ZS
    "Face",       # NSL-ZI
    "Neurocranium", # NA-BR
    "Neurocranium", # NA-FM
    "Face",        # NA-PNS
    "Neurocranium",# BR-PT
    "Neurocranium",# BR-APET
    "Neurocranium",# PT-FM
    "Neurocranium",# PT-APET
    "Neurocranium",# PT-BA
    "Neurocranium",# PT-EAM
    "Face",        # PT-ZYGO
    "Both", # PT-TSP
    "Neurocranium",# FM-ZS
    "Neurocranium",# FM-MT
    "Face",        # ZS-ZI
    "Face",        # ZI-MT
    "Face",        # ZI-ZYGO
    "Face",        # ZI-TSP
    "Face",        # MT-PNS
    "Neurocranium",# PNS-APET
    "Neurocranium",# APET-BA
    "Neurocranium",# APET-TS
    "Neurocranium",# BA-EAM
    "Face",        # EAM-ZYGO
    "Face",        # ZYGO-TSP
    "Neurocranium",# LD-AS
    "Neurocranium",# BR-LD
    "Neurocranium",# OPI-LD
    "Neurocranium",# PT-AS
    "Neurocranium",# JP-AS
    "Neurocranium" # BA-OPI
  )
)

face_traits <- measure_info$Measure[
  measure_info$Module %in% c("Face", "Both")
]

neuro_traits <- measure_info$Measure[
  measure_info$Module %in% c("Neurocranium", "Both")
]

colnames(traits) <- distn

calculate_P_matrices <- function(info, traits, trait_sets) {
  
  genera <- unique(info$GEN)
  genera <- genera[!is.na(genera)]
  
  P_matrices <- list()
  
  for (set_name in names(trait_sets)) {
    
    trait_names <- trait_sets[[set_name]]
    
    P_matrices[[set_name]] <- list()
    
    for (genus in genera) {
      
      rows <- which(info$GEN == genus)
      
      info_genus <- info[rows, ]
      traits_genus <- traits[rows, trait_names, drop = FALSE]
      
      # Combine data
      data_genus <- data.frame(
        traits_genus,
        SPE = info_genus$SPE,
        SEX = info_genus$SEX,
        check.names = FALSE
      )
      
      # Remove individuals with unknown species or sex
      data_genus <- data_genus[
        !is.na(data_genus$SPE) &
          !is.na(data_genus$SEX),
        ,
        drop = FALSE
      ]
      
      # Convert predictors to factors
      data_genus$SPE <- factor(data_genus$SPE)
      data_genus$SEX <- factor(data_genus$SEX)
      
      # Basic requirement: at least two individuals
      if (nrow(data_genus) < 2) {
        warning(
          genus, " (", set_name,
          "): too few individuals (N = ",
          nrow(data_genus), "). Skipping."
        )
        next
      }
      
      # --------------------------------------------------
      # Choose model according to number of species
      # --------------------------------------------------
      
      if (nlevels(data_genus$SPE) == 1) {
        
        # Monospecific genus
        model_type <- "SEX"
        
        model <- try(
          lm(
            as.matrix(
              log(data_genus[, trait_names, drop = FALSE])
            ) ~ SEX,
            data = data_genus
          ),
          silent = TRUE
        )
        
      } else {
        
        # Polyspecific genus
        model_type <- "SPE * SEX"
        
        model <- try(
          lm(
            as.matrix(
              log(data_genus[, trait_names, drop = FALSE])
            ) ~ SPE * SEX,
            data = data_genus
          ),
          silent = TRUE
        )
      }
      
      # Check whether model worked
      if (inherits(model, "try-error")) {
        warning(
          genus, " (", set_name,
          "): model failed. Skipping."
        )
        next
      }
      
      # Calculate P matrix
      P <- try(
        CalculateMatrix(model),
        silent = TRUE
      )
      
      # Check whether CalculateMatrix worked
      if (inherits(P, "try-error")) {
        warning(
          genus, " (", set_name,
          "): CalculateMatrix failed. Skipping."
        )
        next
      }
      
      # Store matrix
      P_matrices[[set_name]][[genus]] <- P
      
      # Report
      cat(
        "Finished:", genus,
        "| Set =", set_name,
        "| Model =", model_type,
        "| N =", nrow(data_genus),
        "| species =", nlevels(data_genus$SPE),
        "| traits =", length(trait_names),
        "\n"
      )
    }
  }
  
  return(P_matrices)
}

trait_sets <- list(
  face = face_traits,
  neuro = neuro_traits
)

P_matrices <- calculate_P_matrices(
  info = info,
  traits = traits,
  trait_sets = trait_sets
)

check_matrix <- function(P) {
  
  P <- as.matrix(P)
  
  # Eigenvalues
  ev <- eigen(P, symmetric = TRUE)$values
  
  # Symmetry
  symmetric <- isTRUE(
    all.equal(P, t(P), tolerance = 1e-8)
  )
  
  # Positive definite
  positive_definite <- all(ev > 0)
  
  # Cholesky check
  chol_ok <- tryCatch(
    {
      chol(P)
      TRUE
    },
    error = function(e) FALSE
  )
  
  # Statistics
  min_ev <- min(ev)
  max_ev <- max(ev)
  n_pos <- sum(ev > 0)
  
  kappa <- if (all(ev > 0)) {
    max_ev / min_ev
  } else {
    NA
  }
  
  logdet <- if (all(ev > 0)) {
    sum(log(ev))
  } else {
    NA
  }
  
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

check_face <- lapply(P_matrices$face, check_matrix)
check_neuro <- lapply(P_matrices$neuro, check_matrix)

check_face_df <- do.call(rbind, lapply(check_face, as.data.frame))
check_neuro_df <- do.call(rbind, lapply(check_neuro, as.data.frame))

#
problematic <- c(
  "Daubentonia",
  "Mirza",
  "Phaner",
  "Arctocebus",
  "Euoticus",
  "Loris",
  "Carlito",
  "Cephalopachus",
  "Tarsius",
  "Prolemur",
  "Rhinopithecus",
  "Allenopithecus"
)

P_matrices$face <- P_matrices$face[
  !names(P_matrices$face) %in% problematic
]

P_matrices$neuro <- P_matrices$neuro[
  !names(P_matrices$neuro) %in% problematic
]

length(P_matrices$face)
length(P_matrices$neuro)

save(
  P_matrices,
  file = "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/P_matrices_face_neuro.RData"
)
