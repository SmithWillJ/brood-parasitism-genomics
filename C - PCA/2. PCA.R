setwd("/file/path")

library(ggplot2)
library(tidyverse)

## Read in the PCA data

pca_data <- read.table("PCA_PLINKPRUNED.eigenvec", header = FALSE)
pca_eigenvals <- read.table("PCA_PLINKPRUNED.eigenval")
location_data <- read.table("GeneFlow_popmap", header = FALSE)
location_column <- location_data[, 2]
pca_data <- cbind(pca_data, Location = location_column) %>%
  rename(Selected = V1)

### Plot PC1 vs PC2

sampleinfo <- read_delim("/file/path/to/sample_info_GeneFlow.csv", delim = ';')

LatitudeInds <- pca_data %>% 
  left_join(sampleinfo %>% mutate(Selected = gsub('_', '', Selected)) %>% dplyr::select(Selected, lat), by = 'Selected') %>%
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
  xlab("PC1 (2.12%)") +
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

#### Plot an identical PCA but with colours of points representing parasitism status.

ParasitePCA_Data <- sampleinfo %>% dplyr::select(Selected, parasit_yes) %>% 
  mutate(parasit_yes = as.character(parasit_yes)) %>%
  mutate(Selected = gsub('_', '', Selected)) %>%
  right_join(pca_data_plot, by = 'Selected')

pca_plot_Parasit <- ggplot(ParasitePCA_Data, aes(V3, V4, colour = as.character(parasit_yes))) +
  geom_point(size = 5) +
  ggtitle("PC1 vs PC2") +
  xlab("PC1 (2.12%)") +
  ylab("PC2 (1.70%)") + scale_colour_manual(values = c("blue","red", "grey"))

pca_plot_Parasit

pca_plot_Parasit + geom_text(aes(label = parasit_yes))

Parasite<- pca_plot_Parasit + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, size = 2),
        axis.title = element_text(size=20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.position = 'none')
        
ggsave('PCA_plot.png', plot = Parasite, path = '~/Desktop', dpi = 600, width = 8, height = 10, unit = 'in')

Parasite

### Plot a map of sites over Europe:

library(sf)
library(dplyr)
library(tidyverse)

map_info <- read_csv('~/file/path/to/map_info.csv') %>%
  dplyr::select(SITE, LONG, LAT, PARASIT, PARASIT_YES) %>%
  st_as_sf(coords = c('LONG', 'LAT'), crs = st_crs(4326)) %>%
  mutate(PARASIT_YES = ifelse(PARASIT_YES == 'Y?', NA, PARASIT_YES)) %>%
  mutate(parasite = ifelse(PARASIT == 0, 'No', 'Yes'))

europe <- st_read('~/file/path/to/ne_50m_land/ne_50m_land.shp') %>%
  st_crop(xmin = -10, ymin = 33.151691,
          xmax = 40, ymax = 63) 

Obj <- read_sf('/file/path/to/RW_map.gpkg') %>%
  st_crop(xmin = -10, ymin = 33.151691,
          xmax = 40, ymax = 63)

Z <- ggplot() +
  geom_sf(data = europe) +
  geom_sf(data=Obj, alpha = 0.2, fill = 'lightgreen') +
  geom_sf(data = map_info, aes(col = parasite), size = 5) +
  scale_color_manual(values = c("blue","red","black"),
                     labels = c("Unparasitised","Parasitised","Unknown")) + theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text = element_text(size = 11),
        legend.position = 'top', legend.title = element_blank())

ggsave('map.png', Z, path = '~/Desktop')

