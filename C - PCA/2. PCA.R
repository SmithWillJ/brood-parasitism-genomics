####################
####################
### PCA plotting ###
####################
####################

### Set working directory ###

setwd("/Users/williamsmith/Desktop/")

### Libraries ###

library(ggplot2)
library(tidyverse)

### Read in data ###

pca_data <- read.table("/filepath/GF_auto_0.5.PCA.eigenvec", header = FALSE)
pca_eigenvals <- read.table("filepath/GF_auto_0.5/GF_auto_0.5.PCA.eigenval")
location_data <- read.table("/filepath/GeneFlow_popmap", header = FALSE)
colnames(location_data) <- c("Selected", "Location")

### Remove underscores from pca_data V1 before joining ###

pca_data <- pca_data %>%
  rename(Selected = V1) %>%
  mutate(Selected = gsub('_', '', Selected)) %>%
  left_join(location_data, by = "Selected")

### Check if Location was joined correctly ###

print(head(pca_data[, c("Selected", "Location")]))

### Prep for plotting PC1 vs PC2 ###
sampleinfo <- read_delim("filepath/sample_info_GeneFlow.csv", delim = ';')

### Create location-latitude mapping ###
location_lat <- pca_data %>%
  left_join(sampleinfo %>% mutate(Selected = gsub('_', '', Selected)) %>% dplyr::select(Selected, lat), by = 'Selected') %>%
  group_by(Location) %>%
  summarise(mean_lat = mean(lat, na.rm = TRUE)) %>%
  arrange(desc(mean_lat))

print(location_lat)

lat_ordered_locations <- location_lat$Location

### Manually set colours for points ###

colour_map <- c(
  "FINW" = "lightblue",
  "FINE" = "darkblue",
  "FINM" = "blue",
  "NO" = "darkgrey",
  "EE" = "violet",
  "DK" = "purple",
  "SE" = "black",
  "LT" = "pink",
  "PL" = "darkgreen",
  "DEDG" = "coral",
  "CZOS" = "burlywood",
  "DEMH" = "darkorange2",
  "CZLU" = "chocolate",
  "SK" = "green",
  "FR" = "brown",
  "ROLS" = "lightgreen",
  "RODD" = "aquamarine",
  "HR" = "brown1",
  "IT" = "orange",
  "ES" = "salmon",
  "TR" = "red",
  "MT" = "deeppink2"
)

### Get PC % explained ###

pca_eigenvals <- read.table("/filepath/GF_auto_0.5.PCA.eigenval")
total_var <- sum(pca_eigenvals$V1)

pc1_pct <- round((pca_eigenvals$V1[1] / total_var) * 100, 2)
pc2_pct <- round((pca_eigenvals$V1[2] / total_var) * 100, 2)
pc3_pct <- round((pca_eigenvals$V1[3] / total_var) * 100, 2)

# Prepare pca_data_plot ###

pca_data_plot <- pca_data %>%
  left_join(sampleinfo %>% mutate(Selected = gsub('_', '', Selected)) %>% dplyr::select(Selected, lat), by = 'Selected') %>%
  mutate(Location = factor(Location, levels = lat_ordered_locations))

print(head(pca_data_plot$Location))

pca_plot1 <- ggplot(pca_data_plot, aes(V3, V4, colour = Location)) +
  geom_point(size = 5) +
  ggtitle("PC1 vs PC2") +
  xlab("PC1 (7.61%)") +
  ylab("PC2 (5.39%)") +
  scale_colour_manual(values = colour_map, breaks = lat_ordered_locations)

pca_plot1

Q<- pca_plot1 +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 2),
        axis.title = element_text(size=20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) + guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))

Q

### Save plot ###

ggsave('PCA_plot_Sites.png', plot = Q, path = '~/Desktop', dpi = 600, width = 13, height = 8, unit = 'in')

### Repeat but colour points by site parasitism status ###

ParasitePCA_Data <- pca_data_plot %>%
  mutate(Selected_noscore = Selected) %>%
  left_join(sampleinfo %>% mutate(Selected = gsub('_', '', Selected)) %>% dplyr::select(Selected, parasit_yes), 
            by = c('Selected_noscore' = 'Selected')) %>%
  mutate(parasit_yes = as.character(parasit_yes))

pca_plot_Parasit <- ggplot(ParasitePCA_Data, aes(V3, V4, colour = parasit_yes)) +
  geom_point(size = 5) +
  ggtitle("PC1 vs PC2") +
  xlab("PC1 (7.61%)") +
  ylab("PC2 (5.39%)") +
  scale_colour_manual(values = c("0" = "blue", "1" = "red"),
                      na.value = "grey")

Parasite<- pca_plot_Parasit +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 2),
        axis.title = element_text(size=20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.position = 'none')

ggsave('PCA_plot_Parasit.png', plot = Parasite, path = '~/Desktop', dpi = 600, width = 8, height = 10, unit = 'in')

Parasite

### Plot a map of sites over Europe ###

library(sf)
library(dplyr)
library(tidyverse)

map_info <- read_csv('/filepath/map_info.csv') %>%
  dplyr::select(SITE, LONG, LAT, PARASIT, PARASIT_YES) %>%
  st_as_sf(coords = c('LONG', 'LAT'), crs = st_crs(4326)) %>%
  mutate(PARASIT_YES = ifelse(PARASIT_YES == 'Y?', NA, PARASIT_YES)) %>%
  mutate(parasite = ifelse(PARASIT == 0, 'No', 'Yes'))

europe <- st_read('/filepath/ne_50m_land/ne_50m_land.shp') %>%
  st_crop(xmin = -10, ymin = 33.151691,
          xmax = 40, ymax = 63)

Obj <- read_sf('/filepath/Site_Map/RW_map.gpkg') %>%
  st_crop(xmin = -10, ymin = 33.151691,
          xmax = 40, ymax = 63)

Z <- ggplot() +
  geom_sf(data = europe) +
  geom_sf(data=Obj, alpha = 0.2, fill = 'lightgreen') +
  geom_sf(data = map_info, aes(col = parasite), size = 5) +
  scale_colour_manual(values = c("blue","red","black"),
                     labels = c("Unparasitised","Parasitised","Unknown")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        axis.text = element_text(size = 11),
        legend.position = 'top', legend.title = element_blank())

ggsave('Map.png', Z, path = '~/Desktop')
