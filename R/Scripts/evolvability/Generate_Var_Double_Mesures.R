library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(openxlsx)
library(MASS)
library(tidyverse)

###################################################################################
# PLATYRRHINI
###################################################################################
df <- openxlsx::read.xlsx("~/Dropbox/Doc/Data/double_measures_platyrrhini.xlsx")
df <- df[,1:85]

# ----------------------------
# 1. Limpar nomes (tirar $)
# ----------------------------
colnames(df) <- gsub("\\$", "", colnames(df))

# Garantir nomes únicos
colnames(df) <- make.unique(colnames(df))

# ----------------------------
# 2. Definir colunas
# ----------------------------
# selecionar apenas colunas numéricas
numeric_cols <- names(df)[sapply(df, is.numeric)]

# remover a coluna 2X da lista de medidas
measure_cols <- setdiff(numeric_cols, c("2X", "PT-TSP", "PT-TSP.1"))

# LOG TRANSFORMATION (ESSENCIAL)
df[measure_cols] <- log1p(df[measure_cols])

# ----------------------------
# 3. Filtrar apenas lados 1 e 2
# ----------------------------
df_filtrado <- df %>%
  filter(`2X` %in% c(1, 2))

# ----------------------------
# 4. Checar IDs problemáticos
# ----------------------------
check_ids <- df_filtrado %>%
  group_by(SPECIES, ID) %>%
  summarise(
    tem_r1 = any(`2X` == 1),
    tem_r2 = any(`2X` == 2),
    n_linhas = n(),
    .groups = "drop"
  )

ids_problematicos <- check_ids %>%
  filter(!tem_r1 | !tem_r2)

# Ver quem está incompleto
print(ids_problematicos)

# ----------------------------
# 5. Transformar para formato longo
# ----------------------------
df_long <- df_filtrado %>%
  pivot_longer(
    cols = all_of(measure_cols),
    names_to = "variavel_raw",
    values_to = "valor"
  ) %>%
  mutate(
    variavel = str_remove(variavel_raw, "\\.\\d+$")
  ) %>%
  filter(variavel != "PTTSP")

# ----------------------------
# 6. Média por lado (r1 e r2)
# ----------------------------
safe_mean <- function(x) {
  if(all(is.na(x))) return(NA)
  mean(x, na.rm = TRUE)
}

resultados <- df_long %>%
  group_by(GENUS, SPECIES, ID, variavel, `2X`) %>%
  summarise(media = safe_mean(valor), .groups = "drop") %>%
  pivot_wider(
    names_from = `2X`,
    values_from = media,
    names_prefix = "r"
  )

# ----------------------------
# 7. (Opcional) manter só IDs completos
# ----------------------------
ids_validos <- check_ids %>%
  filter(tem_r1 & tem_r2)

erros <- resultados %>%
  semi_join(ids_validos, by = c("SPECIES", "ID"))

# ----------------------------
# 8. (Opcional) calcular erro
# ----------------------------
resultados <- resultados %>%
  mutate(diff = r1 - r2)

variancia_erro <- resultados %>%
  group_by(GENUS, SPECIES, variavel) %>%
  summarise(
    var_erro = var(diff, na.rm = TRUE) / 2,
    n = sum(!is.na(diff)),
    .groups = "drop"
  )

# 20 sp de Platyrrhini
p_variancia_erro <- variancia_erro %>%
  filter(n > 1)

saveRDS(p_variancia_erro, "~/Dropbox/Doc/Data/p_variancia_erro.RDS")

###################################################################################
# CATARRHINI
###################################################################################
# -----------------------------
# 1. Ler arquivos dos museus
# -----------------------------

setwd("~/Dropbox/Doc/Data/museus_medidas/")

museus <- c("AIMZU", "AMNH", "FMNH", "MCZ", "MfN", "MNHN-AC", "MNHN-Z",
            "MRAC", "Naturalis", "NHM", "Powell-Cotton", "RBINS", "USNM")

files <- list.files(pattern = "\\.xls$", full.names = TRUE)

c_df_all <- files %>%
  map(~{
    nome_arquivo <- basename(.x)
    museu <- museus[str_detect(nome_arquivo, fixed(museus))]
    
    read_excel(.x, sheet = "transposta") %>%
      mutate(
        MUSEU = museu,
        across(everything(), as.character)
      )
  }) %>%
  bind_rows()

# -----------------------------
# 2. Limpeza básica taxonômica
# -----------------------------
c_df_all <- c_df_all %>%
  filter(!GENUS %in% c("Genero", "genero"),
         !SPECIES %in% c("sp", "0"))

# correções mínimas essenciais
c_df_all <- c_df_all %>%
  mutate(
    GENUS = str_replace(GENUS, "Cercoceb us", "Cercocebus"),
    GENUS = str_replace(GENUS, "Eryhtrocebus", "Erythrocebus"),
    
    SPECIES = str_replace(SPECIES, "erythrostis", "erythrotis"),
    SPECIES = str_replace(SPECIES, "pogonyas", "pogonias"),
    SPECIES = str_replace(SPECIES, "klossi", "klossii"),
    SPECIES = str_replace(SPECIES, "mulleri", "muelleri"),
    SPECIES = str_replace(SPECIES, "nemestrinus", "nemestrina"),
    SPECIES = str_replace(SPECIES, "sylvana", "sylvanus"),
    SPECIES = str_replace(SPECIES, "tesselatum", "tesselatus"),
    SPECIES = str_replace(SPECIES, "cristata", "cristatus"),
    SPECIES = str_replace(SPECIES, "frontatus", "frontata"),
    SPECIES = str_replace(SPECIES, "johni|johnni", "johnii"),
    SPECIES = str_replace(SPECIES, "obscura", "obscurus"),
    SPECIES = str_replace(SPECIES, "phrayrei", "phayrei"),
    SPECIES = str_replace(SPECIES, "rubicundus", "rubicunda"),
    SPECIES = str_replace(SPECIES, "pennanti|pennantti", "pennantii"),
    SPECIES = str_replace(SPECIES, "roxellanae", "roxellana"),
    SPECIES = str_replace(SPECIES, "comatus", "comata")
  )

# -----------------------------
# 3. Padronização geral
# -----------------------------

c_df_all <- c_df_all %>%
  mutate(
    across(c(GENUS, SPECIES, NUMBER, MUSEU),
           ~str_squish(str_to_upper(.x))),
    NUMBER = str_replace_all(NUMBER, "\\.0$", ""),
    MUSEU = recode(MUSEU, "POWELL-COTTON" = "P-C")
  )

# -----------------------------
# 4. Ler msrs
# -----------------------------

msrs <- read.csv(
  "~/Dropbox/Doc/Data/primates_measures/medidas_catarrhini.csv",
  dec = ",", sep = ","
)

msrs <- msrs %>%
  mutate(
    SEX = case_when(
      SEX %in% c("?female") ~ "FEMALE",
      SEX %in% c("?male") ~ "MALE",
      SEX %in% c("0", "", "sexo") ~ NA_character_,
      TRUE ~ SEX
    )
  ) %>%
  drop_na(SEX)

msrs <- msrs %>%
  mutate(
    across(c(GENUS, SPECIES, NUMBER, MUSEUM),
           ~str_squish(str_to_upper(.x))),
    NUMBER = str_replace_all(NUMBER, "\\.0$", "")
  )

# -----------------------------
# 5. Lista de espécies alvo
# -----------------------------

evolvas <- readRDS("~/Dropbox/Doc/Code/evowm/R/Scripts/evolvability/Evolvability.Rds")

species_df <- evolvas$species[1:49] %>%
  tibble(full = .) %>%
  separate(full, into = c("GENUS", "SPECIES"), sep = "_")

clean_taxa <- function(df){
  df %>%
    mutate(
      GENUS   = str_to_upper(str_trim(GENUS)),
      SPECIES = str_to_upper(str_trim(SPECIES))
    )
}

c_df_all <- clean_taxa(c_df_all)
msrs     <- clean_taxa(msrs)

species_df <- evolvas$species[1:49] %>%
  tibble(full = .) %>%
  separate(full, into = c("GENUS", "SPECIES"), sep = "_") %>%
  mutate(
    GENUS = str_to_upper(GENUS),
    SPECIES = str_to_upper(SPECIES)
  )

c_df_all <- semi_join(
  c_df_all,
  species_df,
  by = c("GENUS", "SPECIES")
)

msrs <- semi_join(
  msrs,
  species_df,
  by = c("GENUS", "SPECIES")
)
# -----------------------------
# 6. Traits
# -----------------------------

traits <- c(
  "ISPM","ISNSL","ISPNS","PMZS","PMZI","PMMT",
  "NSLNA","NSLZS","NSLZI","NABR","NAFM","NAPNS",
  "BRPT","BRAPET","PTFM","PTAPET","PTBA","PTEAM",
  "PTZYGO","FMZS","FMMT","ZSZI","ZIMT",
  "ZIZYGO","ZITSP","MTPNS","PNSAPET","APETBA",
  "APETTS","BAEAM","EAMZYGO","ZYGOTSP","LDAS",
  "BRLD","OPILD","PTAS","JPAS","BAOPI"
)

c_df_all[traits] <- lapply(c_df_all[traits], as.numeric)
msrs[traits] <- lapply(msrs[traits], as.numeric)

# -----------------------------
# 7. Média das réplicas
# -----------------------------

c_df_mean <- c_df_all %>%
  group_by(GENUS, SPECIES, MUSEU, NUMBER) %>%
  summarise(
    across(all_of(traits), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# -----------------------------
# 8. Distância
# -----------------------------

calc_dist <- function(x, y){
  valid <- !(is.na(x) | is.na(y))
  if(sum(valid) < 10) return(Inf)
  sqrt(sum((x[valid] - y[valid])^2))
}

# -----------------------------
# 9. Matching
# -----------------------------

matches <- map_lgl(
  1:nrow(c_df_mean),
  function(i){
    
    ref <- as.numeric(c_df_mean[i, traits])
    
    msrs_sub <- msrs %>%
      filter(
        GENUS == c_df_mean$GENUS[i],
        SPECIES == c_df_mean$SPECIES[i],
        MUSEUM == c_df_mean$MUSEU[i]
      )
    
    if(nrow(msrs_sub) == 0) return(FALSE)
    
    dists <- apply(msrs_sub[, traits], 1,
                   function(x) calc_dist(ref, as.numeric(x)))
    
    min(dists, na.rm = TRUE) < 1.5
  }
)

c_df_mean_clean <- c_df_mean[matches, ]

# -----------------------------
# 10. Dataset final
# -----------------------------

c_df_all_clean <- c_df_all %>%
  semi_join(c_df_mean_clean,
            by = c("GENUS","SPECIES","MUSEU","NUMBER"))

# -----------------------------
# 11. Diagnóstico (opcional)
# -----------------------------

best_dists <- map_dbl(
  1:nrow(c_df_mean),
  function(i){
    
    ref <- as.numeric(c_df_mean[i, traits])
    
    msrs_sub <- msrs %>%
      filter(
        GENUS == c_df_mean$GENUS[i],
        SPECIES == c_df_mean$SPECIES[i],
        MUSEUM == c_df_mean$MUSEU[i]
      )
    
    if(nrow(msrs_sub) == 0) return(Inf)
    
    dists <- apply(msrs_sub[, traits], 1,
                   function(x) calc_dist(ref, as.numeric(x)))
    
    min(dists, na.rm = TRUE)
  }
)

c_df_mean$best_dist <- best_dists

unmatched_inf <- c_df_mean %>%
  filter(is.infinite(best_dist))

#

df <- c_df_all_clean  # ou o nome que você estiver usando agora

df <- df %>%
  setNames(make.unique(names(.)))

df_filtrado <- df %>%
  filter(MEASURE...21 %in% c(1, 2))

df_filtrado <- dplyr::select(df_filtrado, -SUBSPECIES)

meta_cols <- c(
  "GENUS", "SPECIES", "NUMBER",
  "SEX", "LOCAL", "COUNTRY",
  "LATITUDE", "LONGITUDE",
  "MEASURE...10", "DATE",
  "COLLECTOR", "ALTITUDE",
  "OBS1","OBS2","OBS3","OBS4","OBS5","OBS6","OBS7",
  "MEASURE...21", "MUSEU"
)

measure_cols <- setdiff(names(df_filtrado), meta_cols)

df_filtrado[measure_cols] <- lapply(df_filtrado[measure_cols], as.numeric)
df_filtrado[measure_cols] <- log1p(df_filtrado[measure_cols])

df_long <- df_filtrado %>%
  pivot_longer(
    cols = all_of(measure_cols),
    names_to = "variavel",
    values_to = "valor"
  )

safe_mean <- function(x){
  if(all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

resultados <- df_long %>%
  group_by(MUSEU, GENUS, SPECIES, NUMBER, variavel, MEASURE...21) %>%
  summarise(
    media = safe_mean(valor),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = MEASURE...21,
    values_from = media,
    names_prefix = "r"
  )

resultados <- resultados %>%
  mutate(diff = r1 - r2)

variancia_erro <- resultados %>%
  group_by(GENUS, SPECIES, variavel) %>%
  summarise(
    n = sum(!is.na(diff)),
    var_erro = ifelse(n > 1, var(diff, na.rm = TRUE) / 2, NA_real_),
    .groups = "drop"
  ) %>%
  filter(n > 1)

c_variancia_erro <- variancia_erro

X <- c_variancia_erro %>%
  distinct(GENUS, SPECIES) %>%
  arrange(GENUS, SPECIES)

#saveRDS(c_variancia_erro, "~/Dropbox/Doc/Data/c_variancia_erro.RDS")

####################################################################
traits <- c(
  "ISPM","ISNSL","ISPNS","PMZS","PMZI","PMMT",
  "NSLNA","NSLZS","NSLZI","NABR","NAFM","NAPNS",
  "BRPT","BRAPET","PTFM","PTAPET","PTBA","PTEAM",
  "PTZYGO","FMZS","FMMT","ZSZI","ZIMT",
  "ZIZYGO","ZITSP","MTPNS","PNSAPET","APETBA",
  "APETTS","BAEAM","EAMZYGO","ZYGOTSP","LDAS",
  "BRLD","OPILD","PTAS","JPAS","BAOPI"
)

# carregar datasets
data_c <- readRDS("~/Dropbox/Doc/Data/c_variancia_erro.RDS")
data_p <- readRDS("~/Dropbox/Doc/Data/p_variancia_erro.RDS")

# padronizar catarrhini
data_c <- data_c %>%
  mutate(
    GENUS = str_to_title(str_to_lower(GENUS)),
    SPECIES = str_to_lower(SPECIES),
    variavel = str_replace_all(variavel, "-", ""),
    variavel = str_to_upper(variavel)
  )

# padronizar platyrrhini
data_p <- data_p %>%
  mutate(
    GENUS = str_to_title(str_to_lower(GENUS)),
    SPECIES = str_to_lower(SPECIES),
    variavel = str_replace_all(variavel, "-", ""),
    variavel = str_to_upper(variavel)
  )

# juntar datasets
data <- bind_rows(data_c, data_p)

# manter apenas traits desejadas
data <- data %>%
  filter(variavel %in% traits)

# criar ordem correta
data$variavel <- factor(data$variavel, levels = traits)

# identificador espécie
data <- data %>%
  mutate(sp = paste(GENUS, SPECIES))

species_list <- unique(data$sp)

error_samples <- list()

for(sp_i in species_list){
  
  df_sp <- data %>%
    filter(sp == sp_i) %>%
    arrange(variavel)
  
  # checar se tem 38 traits
  if(nrow(df_sp) != length(traits)){
    
    message(sp_i, " possui ", nrow(df_sp), " traits")
    next
  }
  
  vars <- df_sp$var_erro
  
  M <- diag(vars)
  
  mu <- rep(0, length(traits))
  
  sims <- MASS::mvrnorm(
    n = 10000,
    mu = mu,
    Sigma = M
  )
  
  sims <- as.data.frame(sims)
  
  colnames(sims) <- traits 
  
  sims$GENUS <- df_sp$GENUS[1]
  sims$SPECIES <- df_sp$SPECIES[1]
  
  sims <- sims %>%
    dplyr::select(GENUS, SPECIES, all_of(traits))
  
  error_samples[[sp_i]] <- sims
}


# juntar
data_all <- bind_rows(data_c2, data_p2) %>%
  filter(variavel %in% traits)

saveRDS(data_all, "~/Dropbox/Doc/Data/var_error.RDS")
saveRDS(error_samples, "~/Dropbox/Doc/Data/error_samples.RDS")

