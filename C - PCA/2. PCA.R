setwd("/file/path/")

library(ggplot2)

## Read in the PCA data

pca_data <- read.table("PCA_PLINKPRUNED.eigenvec", header = FALSE)
pca_eigenvals <- read.table("PCA_PLINKPRUNED.eigenval")
location_data <- read.table("GeneFlow_popmap", header = FALSE)
location_column <- location_data[, 2]
pca_data <- cbind(pca_data, Location = location_column) %>%
  unite('Selected', V1:V2, sep = '_')

### Plot PC1 vs PC2

sampleinfo <- read_csv("/file/path/to/sample_info4.csv")

LatitudeInds <- pca_data %>% 
  left_join(sampleinfo %>% dplyr::select(Selected, lat), by = 'Selected') %>%
  arrange(lat)

Latcolours <- data.frame(Location = factor(unique(LatitudeInds$Location), levels = unique(LatitudeInds$Location[order(-LatitudeInds$lat)])),
                         colrs = c('deeppink2','red','salmon','orange','brown1','aquamarine','lightgreen',
                                   'brown','green','chocolate','coral','burlywood','darkorange2','darkgreen','pink','black','purple','violet', 'darkgrey','blue','darkblue','lightblue'))

pca_data_plot <- pca_data %>%
  left_join(Latcolours, by = 'Location')
pca_data_plot$Location <- factor(pca_data_plot$Location, levels = levels(Latcolours$Location))
color_map <- setNames(pca_data_plot$colrs, pca_data_plot$Location)

pca_plot1 <- ggplot(pca_data_plot, aes(V3, V4, colour = Location)) +
  geom_point(size = 5) +
  ggtitle("PC1 vs PC2") +
  xlab("PC1 (2.13%)") +
  ylab("PC2 (1.70%)") +
  scale_colour_manual(values = color_map)
pca_plot1

pca_plot1 + geom_text(aes(label = Location))

Q<- pca_plot1 + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", 
                                  fill = NA, size = 2),
      axis.title = element_text(size=20, face = 'bold'),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))
ggsave('PCA_plot.png', plot = Q, path = '~/Desktop', dpi = 600, width = 13, height = 8, unit = 'in')

Q

### Plot PC1 vs PC3

pca_plot2 <- ggplot(pca_data, aes(V3, V5, colour = Location)) + geom_point(size = 2) + ggtitle("PC1 vs PC3") +
  xlab("PC1 (2.13%)") + ylab("PC3 (1.64%)") +
  scale_colour_manual(values = c('coral', 'chocolate', 'darkorange2', 'burlywood', 'violet',
                                 'purple', 'salmon', 'darkblue', 'blue', 'lightblue', 'brown', 'brown1',
                                 'orange', 'pink', 'deeppink2', 'darkgrey', 'darkgreen', 'aquamarine',
                                 'lightgreen', 'black', 'green', 'red'))

pca_plot2

pca_plot2 + geom_text(aes(label = Location))

R <- pca_plot2 + theme(panel.background = element_blank(), panel.grid = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 2),
                  axis.title = element_text(size=20, face = 'bold'),
                  axis.text = element_text(size = 15),
                  legend.text = element_text(size = 15),
                  legend.key.size = unit(2.5,"line"),
                  legend.position = "right", legend.title = element_blank(),
                  legend.key = element_rect(fill = NA),
                  legend.background = element_rect(colour="transparent", fill = "transparent"))

### Plot the two figures together

library(gridExtra)
grid.arrange(Q,R, ncol = 2, widths = c(1,1.5))

### Plot a map of sites over Europe:

library(sf)
library(dplyr)
library(tidyverse)

map_info <- read_csv('~/file/path/to/map_info.csv') %>%
  dplyr::select(SITE, LONG, LAT, PARASIT, PARASIT_YES) %>%
  st_as_sf(coords = c('LONG', 'LAT'), crs = st_crs(4326)) %>%
  mutate(PARASIT_YES = ifelse(PARASIT_YES == 'Y?', NA, PARASIT_YES)) %>%
  mutate(parasite = ifelse(PARASIT == 0, 'No', 'Yes'))

europe <- st_read('~/file/path/to/ne_50m_land.shp') %>%
  st_crop(xmin = -10, ymin = 33.151691,
          xmax = 40, ymax = 63) 

Z <- ggplot() +
  geom_sf(data = europe) +
  geom_sf(data = map_info, aes(col = parasite)) +
  geom_sf_text(data = map_info, aes(label = SITE)) +
  scale_color_manual(values = c("red","blue","black"),
                     labels = c("Unparasitised","Parasitised","Unknown")) + theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text = element_text(size = 11),
        legend.position = 'top', legend.title = element_blank())

Z
