# set work directory
setwd("~/Dropbox/Doc/Code/evowm/R/Scripts/")

# load packages
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')

# read all the species files
setwd("~/Dropbox/Doc/Data/wos_mating_systems/Catarrhini/")
temp_c = list.files(pattern="*.xls")
excells_c = lapply(temp_c, read_xls)
names(excells_c)  = gsub(".xls", replacement = "", temp_c)


setwd("~/Dropbox/Doc/Data/wos_mating_systems/Platyrrhini/")
temp_p = list.files(pattern="*.xls")
excells_p = lapply(temp_p, read_xls)
names(excells_p)  = gsub(".xls", replacement = "", temp_p)

temp = c(temp_c, temp_p)

# genus names
gen = unique(sapply(strsplit(temp, "_"), `[`, 1))

# Catarrhini
# define empty vectors
total_c = vector() # total number of articles from WoS
choice_c = vector() # picked articles for mating system research
for(i in 1:length(temp_c)){
  total_c[i] = length(pull(excells_c[[i]][1])) # row number for each species file
  choice_c[i] = sum(excells_c[[i]][,1]) # sum the first column for each species file
}

# Platyrrhini
total_p = vector() # total number of articles from WoS
choice_p = vector() # picked articles for mating system research
for(i in 1:length(temp_p)){
  total_p[i] = length(pull(excells_p[[i]][1])) # row number for each species file
  choice_p[i] = sum(excells_p[[i]][,1]) # sum the first column for each species file
}

# tables of counted papers
all_papers_c = data.frame(names(excells_c), total_c, choice_c, rep("Catarrhini", length(names(excells_c))))
colnames(all_papers_c) = c("Espécies", "Total", "Escolhidos", "Parvorder")

all_papers_p = data.frame(names(excells_p), total_p, choice_p, rep("Platyrrhini", length(names(excells_p))))
colnames(all_papers_p) = c("Espécies", "Total", "Escolhidos", "Parvorder")

all_papers = rbind(all_papers_c, all_papers_p)

write.csv(all_papers, file = "~/Dropbox/Doc/Code/evowm/R/Outputs/Count_Articles.csv")

# basic plots showing the total and picked number distribution per species
hist_total = ggplot(data = all_papers) +
  geom_histogram(aes(Total), alpha = 0.4, color = "black") +
  labs(title = "Frequency of the number of articles found per species in the WoS search", 
       x = "Species", y = "Frequency") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust=0.5))

plot_total = ggplot(data = all_papers) +
  geom_point(aes(x = seq(1,nrow(all_papers),1), y = Total), alpha = 0.4, color = "black", size = 2.5) +
  labs(title = "Total number of articles found per species in the WoS search", 
       x = "Species", y = "Number of Articles") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust=0.5))

##
# Number of articles separating by genus
total_genus = vector()
choice_genus = vector()

for(i in 1:length(gen)){
  index = str_detect(all_papers$Espécies, gen[i])
  for_genus = all_papers[index,]
  
  total_genus[i] = sum(for_genus$Total)
  choice_genus[i] = sum(for_genus$Escolhidos)
}

#  data frame for each genus
all_papers_genus = data.frame(gen, total_genus, choice_genus)
colnames(all_papers_genus) = c("Gênero", "Total", "Escolhidos")

Totais = data.frame(gen, total_genus)
colnames(Totais) = c("gen", "Totals")
Escolhido = data.frame(gen, choice_genus)
colnames(Escolhido) = c("gen", "Choosed")

# Cria o plot e guarda em um objeto
p <- inner_join(Totais, Escolhido) %>% 
  pivot_longer(cols = -gen, names_to = 'Year', names_prefix = 'Count_') %>% 
  ggplot(aes(x = gen, y = value, fill = Year, width = .8)) + 
  geom_col(position = 'dodge') +
  theme_classic(base_size = 18) +
  theme(plot.title = element_text(hjust=0.6)) +
  theme(axis.text.x = element_text(size = 15)) +
  labs(y = "Number of articles", x =  "Primate Genus") +
  scale_fill_discrete(name = "Type") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(" ")

# Salva o plot
ggsave("~/Dropbox/Doc/Code/evowm/R/Outputs/Plot_Articles.png",
       plot = p,          # objeto do plot
       width = 16,        # largura em inches
       height = 8,        # altura em inches
       dpi = 300)         # resolução alta
