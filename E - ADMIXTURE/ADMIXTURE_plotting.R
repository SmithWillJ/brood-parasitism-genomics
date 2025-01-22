setwd("/file/path")
library(pophelper)

remotes::install_github('royfrancis/pophelper')

Admixqlist2 <- readQ("/filepath/to/PLINKPRUNED_EEMS_input_new.2.Q")
Admixqlist3 <- readQ("/filepath/to/PLINKPRUNED_EEMS_input_new.3.Q")
Admixqlist4 <- readQ("/filepath/to/PLINKPRUNED_EEMS_input_new.4.Q")

inds <- read.delim("/filepath/to/GeneFlow_popmap", header = FALSE)

rownames(Admixqlist2[[1]]) <- inds$V1
rownames(Admixqlist3[[1]]) <- inds$V1
rownames(Admixqlist4[[1]]) <- inds$V1

Admixdata2 <- Admixqlist2[[1]] %>% rownames_to_column(var = 'ind') %>%
  pivot_longer(Cluster1:Cluster2, names_to = 'cluster', values_to = 'value') %>%
  mutate(k = 'K = 2')
Admixdata3 <- Admixqlist3[[1]] %>% rownames_to_column(var = 'ind') %>%
  pivot_longer(Cluster1:Cluster3, names_to = 'cluster', values_to = 'value') %>%
  mutate(k = 'K = 3')
Admixdata4 <- Admixqlist4[[1]] %>% rownames_to_column(var = 'ind') %>%
  pivot_longer(Cluster1:Cluster4, names_to = 'cluster', values_to = 'value') %>%
  mutate(k = 'K = 4')

admixdata <- rbind(Admixdata2, Admixdata3, Admixdata4)
popmap <- read_delim('/filepath/to/GeneFlow_popmap', col_names = F) %>%
  rename(ind = X1, pop = X2)
admixdata <- left_join(admixdata, popmap, by = 'ind')


library(ggh4x)
d <- rep('', 188)
dput(d)

sampleinfo <- read_csv("/filepath/to/sample_info4.csv")

LatitudeInds <- admixdata %>% 
  dplyr::select(pop,ind) %>% 
  rename(Selected = ind) %>% 
  left_join(sampleinfo %>% dplyr::select(Selected, lat), by = 'Selected') %>%
  arrange(lat)

Latcolours <- data.frame(pop = factor(unique(LatitudeInds$pop), levels = unique(LatitudeInds$pop[order(-LatitudeInds$lat)])))
admixdata_plot <- admixdata %>%
  left_join(Latcolours, by = 'pop')
admixdata_plot$pop <- factor(admixdata_plot$pop, levels = levels(Latcolours$pop))
admixdata_plot <- admixdata_plot %>%
  mutate(ind = factor(ind, levels = unique(ind)))

(admixture <- ggplot(data = admixdata_plot, aes(x = interaction(ind,pop), y = value, fill = cluster)) + 
    geom_bar(stat = 'identity', position = 'fill') + theme_bw() +
    scale_x_discrete(guide = 'axis_nested') +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_blank(),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = 'bold', hjust = 0)) +
  facet_wrap(~k, ncol = 1))

ggsave('admix.png', plot = admixture, path = '~/Desktop', dpi = 600, width = 13, height = 8, unit = 'in')


plotQMultiline(Admixqlist2, spl = 188, lpp = 1,showindlab = TRUE, outputfilename = "Admix2FINAL",
               useindlab = TRUE, indlabsize = 4.7, 
               dpi = 600, grplabsize = 5, grplabbgcol = "black", grplabcol = "white", indlabspacer = 0.2,
               barsize = 1, imgtype = "jpeg", exportpath = getwd(), grpmean = FALSE,
               clustercol = c("#A765d0","#82AB60"), barbordersize = 0.2, barbordercolour = "white"  ,
               width = 15, height = 5)

plotQMultiline(Admixqlist3, spl = 188, lpp = 1,showindlab = TRUE, outputfilename = "Admix3FINAL",
               useindlab = TRUE, indlabsize = 4.7, 
               dpi = 600, grplabsize = 5, grplabbgcol = "black", grplabcol = "white", indlabspacer = 0.2,
               barsize = 1, imgtype = "jpeg", exportpath = getwd(), grpmean = FALSE,
               clustercol = c("pink","#82AB60","#A765d0"), barbordersize = 0.2, barbordercolour = "white"  ,
               width = 15, height = 5)

plotQMultiline(Admixqlist4, spl = 188, lpp = 1,showindlab = TRUE, outputfilename = "Admix4FINAL",
               useindlab = TRUE, indlabsize = 4.7, 
               dpi = 600, grplabsize = 5, grplabbgcol = "black", grplabcol = "white", indlabspacer = 0.2,
               barsize = 1, imgtype = "jpeg", exportpath = getwd(), grpmean = FALSE,
               clustercol = c("#A765d0","#82AB60","pink","red"), barbordersize = 0.2, barbordercolour = "white"  ,
               width = 15, height = 5)
