# Code to produce analyses and figures for
# Spillover of zoonotic pathogens: A review of reviews 
# by Cecilia Sanchez, Joy Vaz, John Drake
# Code authors: Cecilia Sanchez, Joy Vaz, Eric Marty

# Note that all figures generated in the script were post-processed
# for journal submission by Eric Marty

rm(list = ls())
graphics.off()

# libraries---------------------------------------------------------------------

#install.packages("devtools")
#devtools::install_github('talgalili/dendextend') # dendextend from github

library(tidyverse)
library(magrittr)
library(factoextra)
library(gplots)
library(viridis)
library(data.table)
library(dendextend)
library(reshape2)
library(colorspace)
library(textshape)
library(lsa)
library(wordspace)
library(RColorBrewer)
library(igraph)
library(network)
library(sna)
library(ndtv)
library(circular)
library(cooccur)

# import data ------------------------------------------------------------------
hypo <- read.csv("./Data/Data entry mothersheet - FINAL Hyps.csv", 
                    header = T, as.is = T, encoding = "UTF-8")
tags <- read.csv("./Data/tags.csv", header = T, 
                 as.is = T, encoding = "UTF-8")
colnames(tags) <- c("tag", "theme", "code")

# some merging and cleaning-----------------------------------------------------

# subset hypothesis data
hypo %<>% 
  select(-Timestamp)

# summary of themes
themes <- plyr::count(tags$theme) # Number of tags for each theme
colnames(themes) <- c("theme", "tagcount")
themes$theme <- as.character(themes$theme) # convert from factor to character

# rename hypo columns
names(hypo)[grep("Agri", names(hypo))] <- 
  paste0(themes$theme[grep("AG",themes$theme)],"text")
names(hypo)[grep("Clim", names(hypo))] <- 
  paste0(themes$theme[grep("CL",themes$theme)],"text")
names(hypo)[grep("Huma", names(hypo))] <- 
  paste0(themes$theme[grep("HU",themes$theme)],"text")
names(hypo)[grep("Popu", names(hypo))] <- 
  paste0(themes$theme[grep("PO",themes$theme)],"text")
names(hypo)[grep("Trav", names(hypo))] <- 
  paste0(themes$theme[grep("TR",themes$theme)],"text")
names(hypo)[grep("Soci", names(hypo))] <- 
  paste0(themes$theme[grep("SO",themes$theme)],"text")
names(hypo)[grep("Host", names(hypo))] <- 
  paste0(themes$theme[grep("HO",themes$theme)],"text")
names(hypo)[grep("Path", names(hypo))] <- 
  paste0(themes$theme[grep("PA",themes$theme)],"text")
names(hypo)[grep("Land", names(hypo))] <- 
  paste0(themes$theme[grep("LA",themes$theme)],"text")

# converting hypothesis/tag data to numeric-------------------------------------

# add new columns to populate with info about tags and themes
hypo[, 14:90] <- 0

# name columns for the tags
colnames(hypo)[14:81] <- tags$tag

# name columns for the themes
colnames(hypo)[82:90] <- unique(tags$theme)

# assign a value of 1 if a record includes a tag
# in the corresponding theme column, sum the number of tag occurrences
for(t in 1:nrow(themes)){
  THEME <- as.character(themes$theme[t])
  
  # extract tags from text
  for(i in which(tags$theme == THEME)){
    hypo[, tags$tag[i]][grep(tags$tag[i], hypo[ , paste0(THEME, "text")])] <- 1
  }
  
  # extract themes from tag columns
  hypo[, THEME] <- rowSums(hypo[, with(tags, tag[which(theme == THEME)])])
}

# pull out 0/1 data for tags

codelabs <- tags$code
paperTags <- hypo[, 14:81] 
colnames(paperTags) <- codelabs

# write.csv(hypo, "./Data/FINALcleanedhypo.csv", row.names = FALSE)
# write.csv(paperTags, "./Data/FINALpaperTags.csv", row.names = FALSE)

# Create dend object------------------------------------------------------------

h <- 3
k <- NULL

# Define dendrogram object
set.seed(123) 
dend <- paperTags %>% 
  dist(method = "binary") %>% 
  hclust(method = "ward.D2") %>% 
  as.dendrogram

dend_list <- get_subdendrograms(dend, k = k, h = h)

# palette of colors
clusterpal <- qualitative_hcl(n = length(dend_list), h = c(0, 360), c = 80, 
                              l = 60)

dend <- dend %>% 
  color_branches(col = clusterpal, k = k, h = h)

# Plot dend with clusters colored
par(mfrow = c(1, 1))
plot(dend, main = "Hierarchical Clustering of Hypotheses by Tags", 
     leaflab = "none")

# Add row to hypo to assign cluster numbers to each hyp
clusternum <- dendextend::cutree(dend, k = k, h = h) %>%   
  as.data.frame() %>%  
  setDT(keep.rownames = T) %>% 
  as.tbl() %>% 
  rename('clust' = '.')

hypo_clust <- hypo[, c(1:3, 14:90)] %>% 
  setDT(keep.rownames = T) %>%  
  left_join(clusternum, by = "rn") %>%  
  as.tbl()  

cluster.list <- split(hypo_clust, hypo_clust$clust)

# write .csv files for each cluster --------------------------------------------

# for each cluster, pull out text records from papers and store in .csv
for(i in 1:length(unique(hypo_clust$clust))){
  cluster.i <- as.data.frame(cluster.list[[i]])
  cluster.itext <- cluster.i[, 4]
  write.csv(cluster.itext, paste0("./Data/clusters/clust", i, ".csv"),
            row.names = FALSE)
}

# After reading through hyps in each cluster, assign titles --------------------
cluster.titles <- c("1. Reservoir Hosts",
                    "2. Land Conversion",
                    "3. Global Movement", 
                    "4. Climate & Vectors",
                    "5. Food & Livestock", 
                    "6. Environmental Contact", 
                    "7. Socioeconomics", 
                    "8. Biodiversity & Community",
                    "9. Viral Adaptation")

# re-order to match L to R on dendrogram (Fig. 2)
cluster.titles.reorder <- c("9. Viral Adaptation",
                            "4. Climate & Vectors",
                            "5. Food & Livestock", 
                            "6. Environmental Contact", 
                            "3. Global Movement", 
                            "7. Socioeconomics",
                            "2. Land Conversion",
                            "8. Biodiversity & Community",
                            "1. Reservoir Hosts")

# number of pathways in each cluster (for Table 2)
clust_sizes <- hypo_clust$clust %>%  
  table() %>%  
  as.data.frame() %>% 
  mutate(cluster = cluster.titles)

# Plot dend with clusters colored: FIGURE 2-------------------------------------

par(mfrow = c(1, 1))
plot(dend, main = "Hierarchical Clustering of Hypotheses by Tags", 
     leaflab = "none")

# make colored rectangles around clusters (slow)
fviz_dend(dend, k = k, h = h, k_colors = clusterpal[1:9], rect = T,
          rect_border = clusterpal[1:9], rect_fill = T, lower_rect = 0,
          show_labels = F, lwd = 0.5, sub = "")

# plot subdendrograms-----------------------------------------------------------

# par(mfrow = c(2, 5))
# 
# plot(dend, main = "Original Dendogram", leaflab = "none")
# 
# plotsubdend <- for(i in 1:length(unique(hypo_clust$clust))){
#   cluster.i <- dend_list[[i]]
#   cluster.ititle <- cluster.titles[i]
#   plot(cluster.i, main = cluster.ititle, leaflab = "none")
# }

# multi barplot: FIGURE 3-------------------------------------------------------

resolvedDisp <- read.csv("./Data/FINALpaper_disp_resolved.csv")

names(resolvedDisp) <- c("Paper.ID", "Anthropology", 
                         "Cellular & Molecular Biology", "Ecology & Evolution", 
                         "Environmental Science", "Food & Agriculture", 
                         "Microbiology & Immunology", "Parasitology", 
                         "Public Health & Medicine",
                         "Veterinary Science & Animal Health", "Virology")

disp <- resolvedDisp %>% 
  pivot_longer(`Anthropology`:`Virology`, 
               values_to = "Discipline") %>% 
  filter(Discipline == 1) %>% 
  dplyr::select(-Discipline) %>% 
  rename(PaperDiscipline = name)

hypo_clust %<>% 
  left_join(., disp)

disbyclust <- xtabs(~hypo_clust$PaperDiscipline + hypo_clust$clust)

hyps_per_disp <- hypo_clust %>% 
  group_by(PaperDiscipline) %>% 
  summarise(n = n())

hyps_per_clust <- hypo_clust %>% 
  group_by(clust) %>% 
  summarise(n = n())

for(i in 1:9){
  for(j in 1:10){
    disbyclust[j,i] <- as.numeric(disbyclust[j,i]/hyps_per_clust[i,2]/hyps_per_disp[j,2])
  }
}

frequencies_df <- melt(disbyclust)

frequencies_df %<>%
  rename(Discipline = `hypo_clust$PaperDiscipline`, Cluster = `hypo_clust$clust`,
         Frequency = value) %>%
  mutate_at(vars(Cluster), as.factor) %>% 
  mutate(Cluster = fct_relevel(Cluster, c("9", "4", "5", "6", "3", "7", "2", 
                                          "8", "1")))

dev.off()
ggplot(frequencies_df, aes(x = Cluster, y = Frequency)) +
  geom_bar(stat = "identity", width = 0.4, aes(fill = Cluster)) +
  scale_fill_manual(values = clusterpal) +
  facet_wrap(~Discipline, nrow = 1) +
  coord_flip() + 
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 6, color = "black"),
        axis.text.y = element_text(face = "bold", size = 6, color = "black"),
        axis.title = element_blank(),
        strip.text.x = element_text(face = "bold", size = 6),
        legend.position = "none") +
  scale_x_discrete(labels = c("1" = "Reservoir Hosts", 
                              "8" = "Biodiversity &\nCommunity",
                              "2" = "Land Conversion",
                              "7" = "Socioeconomics",
                              "3" = "Global Movement",
                              "6" = "Environmental Contact",
                              "5" = "Food & Livestock",
                              "4" = "Climate & Vectors",
                              "9" = "Viral Adaptation"))

# heatmap of relative frequency of codes within clusters: FIG 4 (left panel)----

dend_order <- c(9, 4, 5, 6, 3, 7, 2, 8, 1)
dend_order_rev <- rev(dend_order)

cluster.means <- aggregate(paperTags, 
                           by = list(cutree(dend, h = h, k = k)), mean) %>% 
  slice(match(dend_order_rev, Group.1)) %>% 
  column_to_rownames("Group.1")

my_palette <- tritan(sequential_hcl(n = 18, h = 0, c = c(0, NA, NA), 
                                    l = c(10, 98), power = 1.3, rev = FALSE))
gradient_palette <- colorRampPalette(colors = c(my_palette[1:18]), bias = 1) 


dev.off()
heatmap.2(as.matrix(cluster.means),
          col = gradient_palette(20),
          Rowv = T,
          Colv = F,
          RowSideColors = rev(clusterpal),
          labRow = rev(cluster.titles.reorder),
          density.info = "density", 
          denscol = "#00AD9A",
          trace = "none",
          keysize = 0.75,
          key.title = "Key",
          cexCol = 0.94,
          colsep = c(1:67),
          rowsep = c(1:8),
          sepcolor = "grey",
          sepwidth = c(0.05, 0.02),
          margins = c(17, 17))

# heatmap of rel. freq. of codes within disciplines: FIG 4 (right panel)--------

disp.means <- hypo_clust %>%
  dplyr::select(c(agriculture:urbanization, PaperDiscipline)) %>% 
  group_by(PaperDiscipline) %>% 
  summarise_all(mean) %>% 
  select(-PaperDiscipline)

rownames(disp.means) <- c("Anthropology", 
                          "Cellular & Molecular Biology", "Ecology & Evolution", 
                          "Environmental Science", "Food & Agriculture", 
                          "Microbiology & Immunology", "Parasitology", 
                          "Public Health & Medicine",
                          "Veterinary Science & Animal Health", "Virology")

dev.off()
dispx_heatmap <- heatmap.2(as.matrix(disp.means), 
                           col = gradient_palette(20),
                           Rowv = T,
                           Colv = F,
                           #scale = "row",
                           density.info = "density", 
                           denscol = "#00AD9A",
                           trace = "none",
                           keysize = 0.75,
                           key.title = "Key",
                           cexCol = 0.94,
                           colsep = c(1:67),
                           rowsep = c(1:11),
                           sepcolor = "grey",
                           sepwidth = c(0.05, 0.02),
                           margins = c(17, 17))

# network diagram of similarity between disciplines: FIG 5 (right panel)--------

# nodes = paper discipline, edges = cosine similarity

# create links
links <- data.frame(from = rep(paste0("s0", 1:9),
                               times = 9:1),
                    to = c(rep(c(paste0("s0", 2:9), "s10")),
                           rep(c(paste0("s0", 3:9), "s10")),
                           rep(c(paste0("s0", 4:9), "s10")),
                           rep(c(paste0("s0", 5:9), "s10")),
                           rep(c(paste0("s0", 6:9), "s10")),
                           rep(c(paste0("s0", 7:9), "s10")),
                           rep(c(paste0("s0", 8:9), "s10")),
                           c("s09", "s10"),
                           "s10"))
links$weight <- 0

# add resolved PaperDiscipline to hypo
hypo %<>% 
  left_join(., disp)

# create nodes
resolvedDisp2 <- resolvedDisp[, order(colnames(resolvedDisp))]
resolvedDisp2 %<>%
  select(-Paper.ID)

dispSize <- resolvedDisp2 %>% 
  summarise_all(~sum (., na.rm = TRUE)) 

nodes <- data.frame(id = c(paste0("s0", 1:9), "s10"),
                    field = names(resolvedDisp2),
                    size = as.numeric(dispSize[1, ]))

# pull out needed columns: PaperDiscipline and the tags
DTags <- hypo[, c(14:81, 91)]

# sum the number of times a Discipline uses each of the tags
pSums <- DTags %>%
  group_by(PaperDiscipline) %>%
  summarize_all(sum)

# convert tag sums into matrix
pSums2 <- data.matrix(as.data.frame(pSums[, -1]))
# add the paper disciplines back in
row.names(pSums2) <- as.matrix(pSums[, 1])

# calculate cosine similarities between all rows ie disciplines)
cosSim <- dist.matrix(pSums2, method = "cosine", convert = FALSE)

# extract the upper triangular half of the matrix, across by row, no diagonal
links$weight <- t(cosSim)[lower.tri(t(cosSim))]

# helpful link on network diagrams
# http://kateto.net/network-visualization

netPcS <- graph.data.frame(links, nodes, directed = F)

# scales
field_sizePcS <- nodes$size
link_sizePcS <- (links$weight)*10

# layout orientation
lPcS <- layout_in_circle(netPcS)

# function for aligning the labels outside the circle
# https://kieranhealy.org/blog/archives/2011/02/18/aligning-labels-in-circular-igraph-layouts/
radian.rescale <- function(x, start = 0, direction = 1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

lab.locs <- radian.rescale(x = 1:10, direction = -1, start = 0)

# calculating quantiles for coloring the edges
quants <- quantile(links$weight)
colset <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")
E(netPcS)[weight <= quants[2]]$color <- colset[1]
E(netPcS)[weight > quants[2] & weight <= quants[3]]$color <- colset[2]
E(netPcS)[weight > quants[3] & weight <= quants[4]]$color <- colset[3]
E(netPcS)[weight > quants[4]]$color <- colset[4]


dev.off()
par(mar = c(1, 1, 1, 1))

# main graph
plot(netPcS, 
     
     layout = lPcS,      
     
     vertex.shape = "circle", 
     vertex.size = field_sizePcS,
     vertex.color = "white",
     vertex.frame.color = "black",
     
     vertex.label = V(netPcS)$field,
     vertex.label.dist = 2.5,
     vertex.label.degree = lab.locs,
     vertex.label.font = 1, 
     vertex.label.color = "black",
     vertex.label.cex = 1,
     
     edge.width = link_sizePcS,
     edge.curved = 0.2
)

# Network diagram of discipline co-occurrence: FIG 5 (left panel)---------------

# Load and clean data

paper_disp <- read.csv("data/FINALpaper_disp_resolved.csv")

names(paper_disp) <- c("Paper.ID", "Anthropology", 
                       "Cellular & Molecular Biology", "Ecology & Evolution", 
                       "Environmental Science", "Food & Agriculture", 
                       "Microbiology & Immunology", "Parasitology", 
                       "Public Health & Medicine",
                       "Veterinary Science & Animal Health", "Virology")

paper_disp[is.na(paper_disp)] <- 0

disp_paper <- t(paper_disp[2:11])

x <- as.matrix(disp_paper)
x %*% t(x)
y <- x %*% t(x)

# build a graph from the above matrix
g <- graph.adjacency(y, weighted = T, mode = "undirected")

# remove loops
g <- simplify(g)

# set labels and degrees of vertices
V(g)$label <- V(g)$name

V(g)$degree <- colSums(paper_disp)[-1]

# set seed to make the layout reproducible
set.seed(1024)

layout1 <- layout_in_circle(g)

dev.off()
plot(g, layout=layout1)

V(g)$vertex.size <- field_sizePcS # why won't this work
V(g)$label.color <- "black"
V(g)$frame.color <-"black"
V(g)$color <-"white"
E(g)$color <-"black"
egam <- (E(g)$weight)*1.2
E(g)$width <- egam
E(g)$curved <- 0.2

# plot the graph in layout1
plot(g, layout=layout1)
