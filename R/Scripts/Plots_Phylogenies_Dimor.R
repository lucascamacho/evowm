# Pacotes necessários
library(ape)        # para manipulação de árvores
library(phytools)   # para simulação e plot de traços evolutivos
library(viridis)    # para paleta de cores

# datas
mds = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Haplorrhini_MDS_Matings.RDS")
evolvas = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/Evolvability.RDS")
align = read.csv("~/Dropbox/Doc/Code/evowm/R/Outputs/cor_PCS_dimorphism_extant.csv")
medias = readRDS("~/Dropbox/Doc/Code/evowm/R/Outputs/averages_PCS_autovalues_primates.RDS")

#Load and plot the phylogeny
# read and plot phylo tree
filename = "~/Dropbox/Doc/Data/Primates_Dryad_no_scripts/median_tree.tre.nex"
tree = read.nexus(filename)
species = mds$especies
tree = drop.tip(tree, setdiff(tree$tip.label, species))
tree_c = drop.tip(tree, setdiff(tree$tip.label, align$matings.especies[1:38]))
tree_p = drop.tip(tree, setdiff(tree$tip.label, align$matings.especies[39:62]))


# ============================
# 1. Seus dados
# ============================
# tree_p    -> sua árvore de Platyrrhini
dimor = align$dimor
names(dimor) = align$matings.especies
dimor <- dimor[tree_p$tip.label]
sigma2 <- 0.54      # variância do BM
root_state <- 3.91  # estado ancestral estimado


# ============================
# 2. Reconstruir estados ancestrais sob BM
# ============================
anc <- fastAnc(tree_p, dimor, vars = TRUE, CI = TRUE)

# Combinar estados ancestrais e tips para criar paleta
all_values <- c(anc$ace, dimor)

pal <- colorRampPalette(c("blue", "red"))(100)
cols <- pal[as.numeric(cut(all_values,
                           breaks = 100,
                           include.lowest = TRUE))]

# ============================
# 3. Plotar árvore com ramos coloridos
# ============================

# Função para colorir os ramos de acordo com valores médios entre nós/pontas
edge_vals <- numeric(nrow(tree_p$edge))
for(i in 1:nrow(tree_p$edge)){
  parent <- tree_p$edge[i,1]
  child  <- tree_p$edge[i,2]
  
  # Valor do nó ou tip
  val_parent <- if(parent <= length(tree_p$tip.label)) dimor[parent] else anc$ace[parent - length(tree_p$tip.label)]
  val_child  <- if(child <= length(tree_p$tip.label)) dimor[child] else anc$ace[child - length(tree_p$tip.label)]
  
  # Média para colorir o ramo
  edge_vals[i] <- mean(c(val_parent, val_child))
}

# Normalizar e mapear para paleta
edge_cols <- pal[as.numeric(cut(edge_vals,
                                breaks = 100,
                                include.lowest = TRUE))]

#png(
#  filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/tree_dimorfismo_platyrrhini.png",
#  width = 3400,
#  height = 3200,
#  res = 300
#)
pdf(
  file = "~/Dropbox/Doc/Code/evowm/R/Outputs/tree_dimorfismo_platyrrhini.pdf",
  width = 12,   # largura em inches
  height = 10    # altura em inches
)
par(mar = c(5, 5, 4, 8))  # margem direita maior
# Plot
plot(tree_p,
     show.tip.label = TRUE,
     edge.color = edge_cols,
     cex = 1.1,          # ↑ tamanho do nome das espécies
     label.offset = 0.9,# ↑ distância entre o tip e o nome
     lwd = 3             # ↑ grossura dos ramos
)

# Colorir tips com círculos
tipcols <- cols[(length(anc$ace)+1):length(all_values)]
tiplabels(pch = 21, bg = tipcols, cex = 1.8)   # círculos maiores nas pontas
nodelabels(pch = 21, bg = cols[1:length(anc$ace)], cex = 1.4)


val_range <- range(all_values, na.rm = TRUE)

# Legenda
#legend("topright", legend = c("low", "high"), fill = c(pal[1], pal[100]),
#       title = "Dimorfismo Sexual", bty = "n")

add.color.bar(
  leg = 0.4 * max(nodeHeights(tree_p)),  # altura da barra
  cols = pal,
  title = " ",
  lims = val_range,
  digits = 2,
  prompt = FALSE,
  x = max(nodeHeights(tree_p)) * 0.04,
  y = max(nodeHeights(tree_p)) * 0.04
)
dev.off()


#####
dimor = align$dimor
names(dimor) = align$matings.especies
dimor <- dimor[tree_c$tip.label]

integ <- evolvas$Integration
names(integ) <- align$matings.especies
integ <- integ[tree_c$tip.label]


alpha  <- 0.05195238
sigma2 <- 8.38232
beta0  <- 2.481361
beta1  <- 156.0154

anc <- fastAnc(tree_c, dimor, vars = TRUE, CI = TRUE)

# Combinar estados ancestrais e tips para criar paleta
all_values <- c(anc$ace, dimor)

pal <- colorRampPalette(c("blue", "red"))(100)
cols <- pal[as.numeric(cut(all_values,
                           breaks = 100,
                           include.lowest = TRUE))]

# ============================
# 3. Plotar árvore com ramos coloridos
# ============================

# Função para colorir os ramos de acordo com valores médios entre nós/pontas
edge_vals <- numeric(nrow(tree_c$edge))
for(i in 1:nrow(tree_c$edge)){
  parent <- tree_c$edge[i,1]
  child  <- tree_c$edge[i,2]
  
  # Valor do nó ou tip
  val_parent <- if(parent <= length(tree_c$tip.label)) dimor[parent] else anc$ace[parent - length(tree_c$tip.label)]
  val_child  <- if(child <= length(tree_c$tip.label)) dimor[child] else anc$ace[child - length(tree_c$tip.label)]
  
  # Média para colorir o ramo
  edge_vals[i] <- mean(c(val_parent, val_child))
}

# Normalizar e mapear para paleta
edge_cols <- pal[as.numeric(cut(edge_vals,
                                breaks = 100,
                                include.lowest = TRUE))]


#pdf(
#  filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/tree_dimorfismo_catarrhini.pdf",
#  width = 3400,
#  height = 3200,
#  res = 300
#)
pdf(
  file = "~/Dropbox/Doc/Code/evowm/R/Outputs/tree_dimorfismo_catarrhini.pdf",
  width = 12,   # largura em inches
  height = 8    # altura em inches
)
par(mar = c(5, 5, 4, 8))  # margem direita maior
# Plot
plot(tree_c,
     show.tip.label = TRUE,
     edge.color = edge_cols,
     cex = 1.1,          # ↑ tamanho do nome das espécies
     label.offset = 0.9,# ↑ distância entre o tip e o nome
     lwd = 3             # ↑ grossura dos ramos
)

# Colorir tips com círculos
tipcols <- cols[(length(anc$ace)+1):length(all_values)]
tiplabels(pch = 21, bg = tipcols, cex = 1.8)   # círculos maiores nas pontas
nodelabels(pch = 21, bg = cols[1:length(anc$ace)], cex = 1.4)


val_range <- range(all_values, na.rm = TRUE)

# Legenda
#legend("topright", legend = c("low", "high"), fill = c(pal[1], pal[100]),
#       title = "Dimorfismo Sexual", bty = "n")

add.color.bar(
  leg = 0.4 * max(nodeHeights(tree_c)),  # altura da barra
  cols = pal,
  title = " ",
  lims = val_range,
  digits = 2,
  prompt = FALSE,
  x = max(nodeHeights(tree_c)) * 0.01,
  y = max(nodeHeights(tree_c)) * 0.04
)
dev.off()

# SADIQ SULTANSØNN HASSAN

## Caminho do PNG final
png(
  filename = "~/Dropbox/Doc/Code/evowm/R/Outputs/trees_dimorfismo_platy_catarrhini.png",
  width = 5200,   # largura em pixels (proporção 24:10)
  height = 3000,  # altura em pixels
  res = 300       # resolução alta
)

# Layout 1 linha, 2 colunas
par(mfrow = c(1,2), mar = c(5,5,4,8))  # margens iguais

# ============================
# 1. Plot Platyrrhini
# ============================

# Preparar dados
dimor_p <- align$dimor
names(dimor_p) <- align$matings.especies
dimor_p <- dimor_p[tree_p$tip.label]

anc_p <- fastAnc(tree_p, dimor_p, vars = TRUE, CI = TRUE)

all_values_p <- c(anc_p$ace, dimor_p)
pal_p <- colorRampPalette(c("blue", "red"))(100)

# cores
cols_p <- pal_p[as.numeric(cut(all_values_p, breaks = 100, include.lowest = TRUE))]

edge_vals_p <- numeric(nrow(tree_p$edge))
for(i in 1:nrow(tree_p$edge)){
  parent <- tree_p$edge[i,1]
  child  <- tree_p$edge[i,2]
  val_parent <- if(parent <= length(tree_p$tip.label)) dimor_p[parent] else anc_p$ace[parent - length(tree_p$tip.label)]
  val_child  <- if(child <= length(tree_p$tip.label)) dimor_p[child] else anc_p$ace[child - length(tree_p$tip.label)]
  edge_vals_p[i] <- mean(c(val_parent, val_child))
}
edge_cols_p <- pal_p[as.numeric(cut(edge_vals_p, breaks = 100, include.lowest = TRUE))]

# Plot
plot(tree_p,
     show.tip.label = TRUE,
     edge.color = edge_cols_p,
     cex = 1.1,
     label.offset = 0.9,
     lwd = 3,
     main = "Platyrrhini",
     cex.main = 1.8,
     col.main = '#7570b3'
)
tiplabels(pch = 21, bg = cols_p[(length(anc_p$ace)+1):length(all_values_p)], cex = 1.8)
nodelabels(pch = 21, bg = cols_p[1:length(anc_p$ace)], cex = 1.4)

add.color.bar(
  leg = 0.4 * max(nodeHeights(tree_p)),
  cols = pal_p,
  lims = range(all_values_p, na.rm = TRUE),
  digits = 2,
  title = " ",
  cex = 1,
  prompt = FALSE,
  x = max(nodeHeights(tree_p)) * 0.04,
  y = max(nodeHeights(tree_p)) * 0.04
)

# ============================
# 2. Plot Catarrhini
# ============================

dimor_c <- align$dimor
names(dimor_c) <- align$matings.especies
dimor_c <- dimor_c[tree_c$tip.label]

anc_c <- fastAnc(tree_c, dimor_c, vars = TRUE, CI = TRUE)

all_values_c <- c(anc_c$ace, dimor_c)
pal_c <- colorRampPalette(c("blue", "red"))(100)

cols_c <- pal_c[as.numeric(cut(all_values_c, breaks = 100, include.lowest = TRUE))]

edge_vals_c <- numeric(nrow(tree_c$edge))
for(i in 1:nrow(tree_c$edge)){
  parent <- tree_c$edge[i,1]
  child  <- tree_c$edge[i,2]
  val_parent <- if(parent <= length(tree_c$tip.label)) dimor_c[parent] else anc_c$ace[parent - length(tree_c$tip.label)]
  val_child  <- if(child <= length(tree_c$tip.label)) dimor_c[child] else anc_c$ace[child - length(tree_c$tip.label)]
  edge_vals_c[i] <- mean(c(val_parent, val_child))
}
edge_cols_c <- pal_c[as.numeric(cut(edge_vals_c, breaks = 100, include.lowest = TRUE))]

# Plot
plot(tree_c,
     show.tip.label = TRUE,
     edge.color = edge_cols_c,
     cex = 1.1,
     label.offset = 0.9,
     lwd = 3,
     main = "Catarrhini",
     cex.main = 1.8,
     col.main = '#1b9e77'
)
tiplabels(pch = 21, bg = cols_c[(length(anc_c$ace)+1):length(all_values_c)], cex = 1.8)
nodelabels(pch = 21, bg = cols_c[1:length(anc_c$ace)], cex = 1.4)

add.color.bar(
  leg = 0.4 * max(nodeHeights(tree_c)),
  cols = pal_c,
  lims = range(all_values_c, na.rm = TRUE),
  digits = 2,
  title = " ",
  cex = 1,
  prompt = FALSE,
  x = max(nodeHeights(tree_c)) * 0.01,
  y = max(nodeHeights(tree_c)) * 0.04
)

# Fechar arquivo
dev.off()
