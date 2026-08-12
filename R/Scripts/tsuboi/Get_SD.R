setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/")

# ============================================================
# Load data
# ============================================================

load("~/Desktop/Primaset_RawData.RData")
load("~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/vcv.RData")


# ============================================================
# Create lists
# ============================================================

sexual_dimorphism <- list()
genus_means <- list()


# ============================================================
# Loop through genera in VCV
# ============================================================

for (genus in names(vcv)) {
  
  
  # ==========================================================
  # SEXUAL DIMORPHISM
  # ==========================================================
  
  # Individuals from this genus with known sex
  idx <- which(
    info$GEN == genus &
      info$SEX %in% c("M", "F")
  )
  
  
  # Separate sexes
  male_idx <- idx[info$SEX[idx] == "M"]
  female_idx <- idx[info$SEX[idx] == "F"]
  
  
  # Need both sexes for sexual dimorphism
  if (length(male_idx) > 0 && length(female_idx) > 0) {
    
    # Skull measurements
    male_traits <- traits[male_idx, , drop = FALSE]
    female_traits <- traits[female_idx, , drop = FALSE]
    
    
    # Mean of log-transformed measurements
    male_log_gm <- apply(
      male_traits,
      2,
      function(x) mean(log(x), na.rm = TRUE)
    )
    
    female_log_gm <- apply(
      female_traits,
      2,
      function(x) mean(log(x), na.rm = TRUE)
    )
    
    
    # Sexual dimorphism
    # log(Male GM) - log(Female GM)
    SD <- male_log_gm - female_log_gm
    
    
    # Save sexual dimorphism
    sexual_dimorphism[[genus]] <- data.frame(
      Genus = genus,
      N_Males = length(male_idx),
      N_Females = length(female_idx),
      t(SD),
      check.names = FALSE
    )
  }
  
  
  # ==========================================================
  # GENUS MEANS
  # ==========================================================
  
  # ALL individuals from this genus
  genus_idx <- which(info$GEN == genus)
  
  
  # Skull measurements for all individuals
  genus_traits <- traits[genus_idx, , drop = FALSE]
  
  
  # Mean of log-transformed measurements
  genus_log_means <- apply(
    genus_traits,
    2,
    function(x) mean(log(x), na.rm = TRUE)
  )
  
  
  # Save genus means
  genus_means[[genus]] <- data.frame(
    Genus = genus,
    t(genus_log_means),
    check.names = FALSE
  )
}


# ============================================================
# Combine sexual dimorphism data
# ============================================================

sexual_dimorphism_df <- do.call(
  rbind,
  sexual_dimorphism
)

rownames(sexual_dimorphism_df) <- NULL


# ============================================================
# Combine genus means
# ============================================================

genus_means_df <- do.call(
  rbind,
  genus_means
)

rownames(genus_means_df) <- NULL


# ============================================================
# Save files
# ============================================================

save(
  sexual_dimorphism_df,
  file = "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/SD.RData"
)

save(
  genus_means_df,
  file = "~/Dropbox/Doc/Code/evowm/R/Scripts/tsuboi/Genus_Means.RData"
)