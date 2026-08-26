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
  # Get all individuals from this genus
  # ==========================================================
  
  genus_idx <- which(info$GEN == genus)
  
  
  # ==========================================================
  # GENUS MEAN
  # ==========================================================
  
  # Species present in this genus
  species <- unique(info$GSP[genus_idx])
  species <- species[!is.na(species)]
  
  
  # Store species means
  species_means <- list()
  
  
  for (sp in species) {
    
    # Individuals belonging to this species
    sp_idx <- genus_idx[info$GSP[genus_idx] == sp]
    
    sp_traits <- traits[sp_idx, , drop = FALSE]
    
    # Species geometric mean
    # (mean(log(x)) = log(geometric mean))
    sp_log_gm <- apply(
      sp_traits,
      2,
      function(x) mean(log(x), na.rm = TRUE)
    )
    
    species_means[[sp]] <- sp_log_gm
  }
  
  
  # Combine species means
  species_means <- do.call(
    rbind,
    species_means
  )
  
  
  # Mean across species
  # Every species gets equal weight
  genus_log_mean <- colMeans(
    species_means,
    na.rm = TRUE
  )
  
  
  genus_means[[genus]] <- data.frame(
    Genus = genus,
    N_Species = nrow(species_means),
    t(genus_log_mean),
    check.names = FALSE
  )
  
  
  # ==========================================================
  # SEXUAL DIMORPHISM
  # ==========================================================
  
  species_SD <- list()
  
  
  for (sp in species) {
    
    # Individuals from this species with known sex
    sp_idx <- genus_idx[
      info$GSP[genus_idx] == sp &
        info$SEX[genus_idx] %in% c("M", "F")
    ]
    
    
    # Separate sexes
    male_idx <- sp_idx[info$SEX[sp_idx] == "M"]
    female_idx <- sp_idx[info$SEX[sp_idx] == "F"]
    
    
    # Need both sexes
    if (length(male_idx) > 0 && length(female_idx) > 0) {
      
      male_traits <- traits[male_idx, , drop = FALSE]
      female_traits <- traits[female_idx, , drop = FALSE]
      
      
      # Species-level male geometric mean
      male_log_gm <- apply(
        male_traits,
        2,
        function(x) mean(log(x), na.rm = TRUE)
      )
      
      
      # Species-level female geometric mean
      female_log_gm <- apply(
        female_traits,
        2,
        function(x) mean(log(x), na.rm = TRUE)
      )
      
      
      # Species-level sexual dimorphism
      SD <- male_log_gm - female_log_gm
      
      
      # Save
      species_SD[[sp]] <- SD
    }
  }
  
  
  # ==========================================================
  # GENUS-LEVEL SEXUAL DIMORPHISM
  # ==========================================================
  
  if (length(species_SD) > 0) {
    
    species_SD <- do.call(
      rbind,
      species_SD
    )
    
    
    # Mean sexual dimorphism across species
    # Every species gets equal weight
    genus_SD <- colMeans(
      species_SD,
      na.rm = TRUE
    )
    
    
    sexual_dimorphism[[genus]] <- data.frame(
      Genus = genus,
      N_Species = nrow(species_SD),
      t(genus_SD),
      check.names = FALSE
    )
  }
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