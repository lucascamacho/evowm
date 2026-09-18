##### LER #####

setwd("~/Downloads/")

dados = read.csv("nhm.csv")

index = which(dados$order == "Primates")
dados = dados[index,]

new = data.frame(dados$sex, dados$lifeStage, dados$verbatimLocality, dados$genus, dados$specificEpithet, 
                 dados$species)

write.table(dados, "NHM_Primates.txt")




### CRIAR PLANILHA
setwd("~/Downloads/")

# Carregar pacotes necessários
library(tidyverse)
library(readxl)
library(writexl)
library(tidyr)

# Carregar os pacotes necessários
if (!require(tidyr)) install.packages("tidyr")
library(tidyr)

dados = read_xlsx("amostras.xlsx")

# Garantir que as colunas "Machos" e "Femeas" estão no formato numérico
dados[, -1] <- lapply(dados[, -1], as.numeric)

# Converter os dados para formato longo
dados_long <- pivot_longer(dados, 
                           cols = -Especies, 
                           names_to = c("Museu", ".value"), 
                           names_sep = "_")

# Conferir o novo formato
head(dados_long)

# Prioridade dos museus
prioridade <- c("RMCA", "AMNH", "NHM", "MFN", "MNHN", "RBINS", "FMNH")

# Reordenar os museus pela prioridade
dados_long <- dados_long[order(match(dados_long$Museu, prioridade)), ]

# Função para calcular a amostragem por espécie
amostrar_especie <- function(df) {
  # Checar se existe um museu com pelo menos 20 machos e 20 fêmeas
  museu_completo <- df[df$Machos >= 20 & df$Femeas >= 20, ]
  
  # Se existir, amostrar apenas esse museu
  if (nrow(museu_completo) > 0) {
    return(museu_completo[1, , drop = FALSE])
  }
  
  # Caso contrário, seguir a estratégia de acumular amostras
  total_machos <- 0
  total_femeas <- 0
  amostragem <- list()
  
  for (i in seq_len(nrow(df))) {
    disponivel_machos <- df$Machos[i]
    disponivel_femeas <- df$Femeas[i]
    
    # Coletar o máximo possível sem exceder 20
    coletar_machos <- min(20 - total_machos, disponivel_machos)
    coletar_femeas <- min(20 - total_femeas, disponivel_femeas)
    
    if (coletar_machos > 0 || coletar_femeas > 0) {
      linha <- df[i, ]
      linha$Machos <- coletar_machos
      linha$Femeas <- coletar_femeas
      amostragem[[length(amostragem) + 1]] <- linha
    }
    
    # Atualiza os totais
    total_machos <- total_machos + coletar_machos
    total_femeas <- total_femeas + coletar_femeas
    
    # Parar se atingir 20 machos e 20 fêmeas
    if (total_machos >= 20 && total_femeas >= 20) break
  }
  
  return(do.call(rbind, amostragem))
}

# Aplicar a função a cada espécie
resultado <- do.call(rbind, lapply(split(dados_long, dados_long$Especies), amostrar_especie))

# Conferir o resultado final
head(resultado)

# Salvar em um novo arquivo Excel
write.csv(resultado, "resultado_amostragem.csv", row.names = FALSE)

##### PEGAR SP VALIDAS #####
library(dplyr)

resumo_especies <- resultado %>%
  group_by(Especies) %>%
  summarise(Machos = sum(Machos), Femeas = sum(Femeas), .groups = "drop")

# Filtrar espécies com no mínimo 20 machos e 20 fêmeas
especies_validas <- resumo_especies %>%
  filter(Machos >= 18 & Femeas >= 18)

# Exibir os nomes das espécies que atendem ao critério
print(especies_validas$Especies)
