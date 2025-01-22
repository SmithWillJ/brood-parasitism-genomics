#########################################################################
#########################################################################
#### (p)RDA ANALYSIS - Reed Warbler parasitism - Will Smith 2024 ########
#########################################################################
### adapted from: https://popgen.nescent.org/2018-03-27_RDA_GEA.html ####
#########################################################################
#########################################################################

### Check BiocManager is installed. If not, install it ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

### Use Biocmanager to install LEA ###

BiocManager::install("LEA")

### Load psych, vegan, and LEA ###

library(tidyverse)
library(psych)
library(vegan)
library(LEA)

### Read in the genomic data. ###

snp<-read.geno("file/path/to/populations.snps.geno")

### Print dimensions of the genomic data ###

dim(snp)

### Read in the environmental data (removing the first column, and printing the structure) ###

env <- read.csv("/file/path/to/sample_info2.csv")
env <- env[,-1]
str(env)

### Change column which is individual names to character type, so it's common to env and genomic datasets ###

env$Selected <- as.character(env$Selected)
sel <- env$Selected
rownames(snp) <- sel

### Check if the row names are the same in the genomic and environmental datasets ###

identical(rownames(snp), env[,1])

### Print the structure of the genomic dataset ###

str(snp)

### Assess pairs panels for the environmental variables ###

pairs.panels(env[,6:24], scale =T)

### Subset the environmental data to keep the variables of interest ###

pred <- subset(env,select = c(parasit, lat, bio4, bio8, bio13, bio15))

### Perform RDA analysis with genomic data as response variable and environmental data as predictors ###

test.rda <- rda(snp ~ ., data=pred, scale=T, set.seed(1))

### Print the RDA analysis results ###

test.rda
RsquareAdj(test.rda)

### This constrained ordination explains a tiny amount of total variation, ###
### But this is expected - most SNPs are going to be neutral ###
### or correlate with other environmental variables ###

### Summarise eigenvalues ###

summary(eigenvals(test.rda,  model = "constrained"))

### Visualise a screeplot ###

screeplot(test.rda)

### Perform an ANOVA on the RDA to test the significance of all constraints ###
### Global RDA is significant, p = 0.001 ###

signif.full <- anova.cca(test.rda, parallel = getOption("mc.cores"))
signif.full

### Perform an ANOVA to check each constrained axis for significance ###

signif.axis <- anova.cca(test.rda, by = "axis", parallel = getOption("mc.cores"))
signif.axis

### Check Variance Inflation Factors of predictors #######
### multicollinearity should not be an issue ###

vif.cca(test.rda)

### We will now plot the RDA ###
### Plot RDA1 & RDA2, scaling SNP/individual scores by ###
### the square root of the eigenvalues. SNPs are grey, and individuals are coloured ###
### The ordination axes are linear combinations of the predictor variables. ###

snp_gg_data <- as.data.frame(test.rda$CCA$v)
ind_gg_data <- as.data.frame(test.rda$CCA$wa) %>%
  rownames_to_column(var = 'name')
loadings <- as.data.frame(test.rda$CCA$biplot)
locations <- env %>% dplyr::select(Selected, lon, lat) %>%
  group_by(lon, lat) %>%
  tally() %>%
  ungroup() %>%
  mutate(location = row_number()) %>%
  mutate(location = str_pad(location, width = 2, pad = 0))
name_loc <- env %>% dplyr::select(Selected, lon, lat) %>% left_join(locations, by = c('lon', 'lat')) %>%
  rename(name = Selected)
ind_gg_data <- left_join(ind_gg_data, name_loc , by = 'name')

PLOT_A <- ggplot() +
  geom_point(data = snp_gg_data, aes(x = RDA1, y = RDA2), col = 'grey', pch = 3) +
  geom_point(data = ind_gg_data, aes(x = RDA1, y = RDA2, col = as.character(location))) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = RDA1*0.1, yend = RDA2*0.1), 
                                    arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  geom_text(data = loadings, aes(x = RDA1*0.13, y = RDA2*0.13), label = c('Parasitism', 'Latitude', 'BIO4', 'BIO8', 'BIO13', 'BIO15'), fontface = 'bold') +
  scale_color_manual(values = c('brown', 'salmon', 'darkorange2', 'burlywood', 'darkgrey',
                                'violet', 'deeppink2', 'orange', 'coral', 'darkgreen', 'brown1','green', 'pink',
                                'lightblue', 'lightgreen', 'purple', 'blue', 'darkblue', 'aquamarine','red'),
                       labels = c("FR", "ES", "DEDG", "DEMH", "NO", "DK", "MT", "IT", "CZLU", "PL", "HR", "SK", "LT", "FINW", "ROLS", "EE", "FINM", "FINE", "RODD","TR")) +
  xlab("RDA1 (23.90%)") + ylab("RDA2 (17.04%)") + theme_bw() +
  theme(legend.position = 'none') + theme(panel.grid = element_blank()) + theme(axis.text = element_text(size = 11),
                                                                                      axis.title = element_text(size = 20))

PLOT_A  

### Assess SNP loadings in ordination space

load.rda <- scores(test.rda, choices = c(1:5), display = 'species')
load.rda

### Assess the loadings on each RDA axis. SNPs located at the tails of these ### 
### distributions will be more likely to be under selection as a function of ###
### these predictors (or correlated predictors!) ###

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
hist(load.rda[,4], main="Loadings on RDA4")
hist(load.rda[,5], main="Loadings on RDA5")

### Define a function to get values in vector 'x' which are outside ###
### of 'z' standard deviations from the mean ###

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)          
  x[x < lims[1] | x > lims[2]]               
}

### Apply this function to the loadings of the first five RDA axes, ### 
### which we know to be significant ###

cand1 <- outliers(load.rda[,1],3)
cand1
cand2 <- outliers(load.rda[,2],3)
cand2
cand3 <- outliers(load.rda[,3],3)
cand3
cand4 <- outliers(load.rda[,4],3)
cand4
cand5 <- outliers(load.rda[,5],3)
cand5

### Get the total number of outliers we have identified ###

ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5)
ncand

### Creates dataframes for each axis with information about the SNPs and the ###
### loadings of potential outliers ###

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))

### Give the dataframes column names ###

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4) <- colnames(cand5) <- c("axis","snp","loading")

### Combine dataframes, and make sure that the 'snp' column is of type 'character' ###

cand <- rbind(cand1, cand2, cand3, cand4, cand5)
cand$snp <- as.character(cand$snp)

### Get the correlations of each candidate SNP with our predictors ###

foo <- matrix(nrow=(ncand), ncol=6)  
colnames(foo) <- c("parasit", "lat", "bio4", "bio8", "bio13", "bio15")
cand$snp <- gsub('col', '', cand$snp)
cand$snp <- as.numeric(cand$snp)
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

### Combine the foo and cand datasets to get our candidate file ###

cand <- cbind.data.frame(cand,foo)  
head(cand)

### Look for duplicated SNPs ###

length(cand$snp[duplicated(cand$snp)])

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1

table(foo[foo[,1]==2,2]) #  no duplicates on axis 2

table(foo[foo[,1]==3,2]) # duplicates on axis 3

table(foo[foo[,1]==4,2]) # duplicates on axis 4

table(foo[foo[,1]==5,2]) # duplicates on axis 5

## Remove duplicates ###

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections to retain only unique SNPs

### Now let's look at which predictors each candidate SNP is most strongly correlated with ###
### This code iterates through the rows of 'cand', identifies the variables with the ###
### maximum correlation, and stores these results in columns 9 and 10 ###

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,10] <- names(which.max(abs(bar[4:8]))) # gives the variable
  cand[i,11] <- max(abs(bar[4:8]))              # gives the correlation
}

colnames(cand)[10] <- "predictor"
colnames(cand)[11] <- "correlation"

### SO HOW MANY CANDIDATES FOR EACH VARIABLE DO WE HAVE? ###

table(cand$predictor)

### 41 are associated with parasitism ###

ParaSNPs <- cand[cand$predictor == "parasit","snp"]
values <- unique(ParaSNPs)
values_para <- unique(ParaSNPs)
values

LatSNPs <- cand[cand$predictor == "lat","snp"]
valuesLat <- unique(LatSNPs)
values_Lat <- unique(LatSNPs)
valuesLat

Bio4SNPs <- cand[cand$predictor == "bio4","snp"]
valuesBio4 <- unique(Bio4SNPs)
values_Bio4 <- unique(Bio4SNPs)
valuesBio4

Bio8SNPs <- cand[cand$predictor == "bio8","snp"]
valuesBio8 <- unique(Bio8SNPs)
values_Bio8 <- unique(Bio8SNPs)
valuesBio8

Bio13SNPs <- cand[cand$predictor == "bio13","snp"]
valuesBio13 <- unique(Bio13SNPs)
values_Bio13 <- unique(Bio13SNPs)
valuesBio13

### We will finish off by looking at the RDA plots and working out where ###
### these SNPs fall out in the ordination space ###

### Get the colours ###

sel <- cand$snp
env <- cand$predictor
env[env=="parasit"] <- 'red'
env[env=="lat"] <- 'grey'
env[env=="bio4"] <- 'forestgreen'
env[env=="bio8"] <- 'green'
env[env=="bio13"] <- 'lightgreen'
env[env=="bio15"] <- 'darkgreen'

### Assign colours to candidate SNPs, and make non-candidate SNPs drab ###

col.pred <- rownames(test.rda$CCA$v)
sel <- paste0('col', sel)

for (i in 1:length(sel)) {           
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("col",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('red','grey','forestgreen','green','lightgreen','darkgreen')

str(bg)
str(env)

### Generate the plot for RDA 1 versus RDA 2 ###

plot(test.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(test.rda, display="species", pch=21, cex=1, col="gray32", bg = col.pred, scaling=3)
points(test.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(test.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("parasit","lat","bio4","bio8","bio13","bio15"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

### Plot for RDA 1 versus RDA 3 ###

plot(test.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(test.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(test.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(test.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("parasit","lat","bio4","bio8","bio13","bio15"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

#########################################################################
#########################################################################
#### PARTIAL RDA ########################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#### First, we will run a PCA to characterise genomic variation in the SNP dataset #####
### Perform a PCA on the genomic data only ###

gen_data <- prcomp(snp, scale. = T)
screeplot(gen_data)
gen_data

#### Seems sensible to retain PC1 as a proxy for population genetic structure #####
PCs <- scores(gen_data, choices = c(1:1), display = "sites", scaling=0)
library(tidyverse)
PCs <- as.data.frame(PCs)
PCs <- PCs%>%rownames_to_column(var = ("Selected"))
env_gen <- left_join(as.data.frame(env),PCs,by=("Selected"))
dim(env_gen)

### Full model ###

pRDAfull <- rda(snp ~ PC1 + parasit + lat + bio4 + bio8 + bio13 + bio15, env_gen)
RsquareAdj(pRDAfull)
anova(pRDAfull)

### Pure population genetic structure model ###

pRDAstruct <- rda(snp ~ PC1 + Condition(parasit + lat + bio4 + bio8 + bio13 + bio15), env_gen)
RsquareAdj(pRDAstruct)
anova(pRDAstruct)

#### Genome scan using pRDA #####

RDA_env <- rda(snp ~ parasit + lat + bio4 + bio8 + bio13 + bio15 + Condition(PC1), env_gen, scale=T, set.seed(1))
screeplot(RDA_env, main="Eigenvalues of contrained axes, pRDA")
RDA_env
RsquareAdj(RDA_env)
summary(eigenvals(RDA_env,  model = "constrained"))
signif.full <- anova.cca(RDA_env, parallel = getOption("mc.cores"))
signif.full
signif.axis <- anova.cca(RDA_env, by = "axis", parallel = getOption("mc.cores"))
signif.axis
vif.cca(RDA_env)
plot(RDA_env, scaling = 3)

#### Repeat plotting as for RDA ####

snp_gg_data <- as.data.frame(RDA_env$CCA$v)
ind_gg_data <- as.data.frame(RDA_env$CCA$wa) %>%
  rownames_to_column(var = 'name')
loadings <- as.data.frame(RDA_env$CCA$biplot)
locations <- env %>% dplyr::select(Selected, lon, lat) %>%
  group_by(lon, lat) %>%
  tally() %>%
  ungroup() %>%
  mutate(location = row_number()) %>%
  mutate(location = str_pad(location, width = 2, pad = 0))
name_loc <- env %>% dplyr::select(Selected, lon, lat) %>% left_join(locations, by = c('lon', 'lat')) %>%
  rename(name = Selected)
ind_gg_data <- left_join(ind_gg_data, name_loc , by = 'name')

PLOT_B <- ggplot() +
  geom_point(data = snp_gg_data, aes(x = RDA1, y = RDA2), col = 'grey', pch = 3) +
  geom_point(data = ind_gg_data, aes(x = RDA1, y = RDA2, col = as.character(location))) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = RDA1*0.1, yend = RDA2*0.1), 
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  geom_text(data = loadings, aes(x = RDA1*0.13, y = RDA2*0.13), label = c('Parasitism', 'Latitude', 'BIO4', 'BIO8', 'BIO13', 'BIO15'), fontface = 'bold') +
  scale_color_manual(values = c('brown', 'salmon', 'darkorange2', 'burlywood', 'darkgrey',
                                'violet', 'deeppink2', 'orange', 'coral', 'darkgreen', 'brown1','green', 'pink',
                                'lightblue', 'lightgreen', 'purple', 'blue', 'darkblue', 'aquamarine','red'),
                     labels = c("FR", "ES", "DEDG", "DEMH", "NO", "DK", "MT", "IT", "CZLU", "PL", "HR", "SK", "LT", "FINW", "ROLS", "EE", "FINM", "FINE", "RODD","TR")) +
  xlab("RDA1 (19.12%)") + ylab("RDA2 (18.09%)") + theme_bw() +
  theme(legend.position='none') + theme(panel.grid = element_blank()) + theme(axis.text = element_text(size = 11),
                                                                                      axis.title = element_text(size = 20))
PLOT_B

#### Identify outliers as for RDA #####

load.rda <- scores(RDA_env, choices = c(1:5), display = 'species')
load.rda
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
hist(load.rda[,4], main="Loadings on RDA4")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3)
cand1
cand2 <- outliers(load.rda[,2],3)
cand2
cand3 <- outliers(load.rda[,3],3)
cand3
cand4 <- outliers(load.rda[,4],3)
cand4
ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4)
ncand
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3, cand4)
cand$snp <- as.character(cand$snp)
foo <- matrix(nrow=(ncand), ncol=6)  # 6 columns for 6 predictors
colnames(foo) <- c("parasit", "lat", "bio4", "bio8", "bio13", "bio15")
cand$snp <- gsub('col', '', cand$snp)
cand$snp <- as.numeric(cand$snp)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snp[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)
length(cand$snp[duplicated(cand$snp)]) ### two duplicates
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #  duplicates on axis 2
table(foo[foo[,1]==3,2]) # duplicates
table(foo[foo[,1]==4,2]) # duplicates

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,10] <- names(which.max(abs(bar[4:8]))) # gives the variable
  cand[i,11] <- max(abs(bar[4:8]))              # gives the correlation
}

colnames(cand)[10] <- "predictor"
colnames(cand)[11] <- "correlation"
table(cand$predictor)

ParaSNPsB <- cand[cand$predictor == "parasit","snp"]
valuesB <- unique(ParaSNPsB)
valuesB
values

all_values <- as.data.frame(c(values, valuesB)) 
write_delim(all_values, '~/Desktop/candidate_snps.txt', delim = ' ')

shared_values <- intersect(values, valuesB)
shared_values
values_only <- setdiff(values,valuesB)
values_only
valuesB_only <-setdiff(valuesB, values)
valuesB_only

cat("Number of values shared between the two lists:", length(shared_values), "\n")
cat("Number found in values only:", length(values_only), "\n")
cat("Number found in valuesB only:", length(valuesB_only), "\n")

########

## Manhattan Plot of LFMM results - with (p)RDA SNPs also highlighted #####

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read in the data
ooglie_booglies <- read.csv('/Users/williamsmith/Desktop/Reed_Warblers_Parasitism/Reed_Warbler_Figures/Genes_and_polygenic_scores/ooglie_booglies.csv')

# Define model and candidates
model <- "parasit"
candidates <- c("2735", "3027", "5216", "5970", "8689", "10153", "11584","12411",
                "12616","12898","2753","6811","9252","10611","11818","363","392",
                "1993", "2095", "2133", "2963","3515", "3546", "3932","4026","4076",
                "4764", "5110","5890","5910","6225", "8538","9665","9769", "9831","11110","11114","11761", "11767","12547","13160","9116","3512","3828","3973",
                "4418", "12649", "9806","4602")
candidates

# Create a function to label selected points
label_points <- function(df, x_col, y_col, label_col, labels) {
  df <- df %>%
    filter(label_col %in% labels)
  df <- df[order(df[[x_col]]), ]
  df$x_pos <- df[[x_col]]
  df$y_pos <- df[[y_col]]
  return(df)
}

library(tidyverse)

# Add positional index to the dataset; then add alternating colour for chromosome;
data <- read_csv('/file/path/to/data.csv') %>%
  separate(V1, into = c('super', 'chrom'), sep = '_') %>%
  mutate(chrom = str_pad(chrom, width = 2, pad = 0)) %>%
  mutate(chrom = ifelse(chrom == '0Z', 30, ifelse(chrom == '0W', 31, chrom))) %>%
  arrange(chrom) %>%
  mutate(index = row_number()) %>%
  filter(super == 'SUPER') %>%
  group_by(chrom) %>%
  mutate(min = min(index), max = max(index)) %>%
  mutate(centre = (max + min)/2) %>%
  mutate(candidate = ifelse(V4 %in% candidates, 'yes', 'no')) %>%
  mutate(magic = ifelse(V4 %in% c("10153", "3027", "8689", "9831", "5890", "3515", "392", "11761", "363", "3546", "11110"), 'magic', 'nie'),
         bayes = ifelse(V4 %in% c("9806","4602"), 'bayes', 'nie'))

# Create the plot
gene_names <- read_delim('/file/path/to/just_genes.txt') %>% rename(V4 = snpID)
ooglie_booglies <- data %>% left_join(gene_names, by = 'V4') 
top_hits <- data %>% filter(candidate == 'yes') %>% ungroup() %>% slice_max(-log10(parasit), n = 10)

library(ggtext)
PLOT_C <- ggplot(data = data, aes(x = index, y = -log10(parasit), col = (as.numeric(chrom) %% 2 == 0))) + 
  geom_point() + labs(y = "-log<sub>10</sub>(<i>p</i>)") +
  geom_point(data = data %>% filter(candidate == 'yes'), color = 'blue') +
  geom_point(data = data %>% filter(candidate == 'yes') %>% ungroup() %>% slice_max(-log10(parasit), n = 11), size = 3, shape = 1, col = 'black') +
  geom_point(data = data %>% filter(magic == 'magic'), color = 'red') +
  geom_point(data = data %>% filter(bayes == 'bayes'), color = 'magenta') +
  ggrepel::geom_text_repel(data = data %>% filter(candidate == 'yes') %>% ungroup() %>% slice_max(-log10(parasit), n = 26), aes(label = Name), col = 'black', size =20/.pt) +
  scale_color_manual(values = c('grey84', 'grey44'), guide = 'none') +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 1, 2, 3, 4, 5, 6)) +
  scale_x_continuous(name = 'Chromosome',label = c(1,2,3,4,5,6,7,8,9,10,11,12,
                                                   13,14,15,16,17,18,19,20,
                                                   21,22,23,'','','','','','','Z', 'W'), breaks = unique(data$centre)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_markdown(), axis.text = element_text(size = 11), 
        axis.title = element_text(size = 20)) 

PLOT_C

### LEGEND PLOT
PLOT_LEG <- ggplot() +
  geom_point(data = snp_gg_data, aes(x = RDA1, y = RDA2), col = 'grey', pch = 3) +
  geom_point(data = ind_gg_data, aes(x = RDA1, y = RDA2, col = as.character(location))) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = RDA1*0.1, yend = RDA2*0.1), 
               arrow = arrow(length = unit(1/2, "picas")), color = "black") +
  geom_text(data = loadings, aes(x = RDA1*0.13, y = RDA2*0.13), label = c('Parasitism', 'BIO1', 'BIO4', 'BIO8', 'BIO13', 'BIO17'), fontface = 'bold') +
  scale_color_manual(name = 'Location',values = c('brown', 'salmon', 'darkorange2', 'burlywood', 'darkgrey',
                                                  'violet', 'deeppink2', 'orange', 'coral', 'darkgreen', 'brown1','green', 'pink',
                                                  'lightblue', 'lightgreen', 'purple', 'blue', 'darkblue', 'aquamarine','red'),
                     labels = c("FR", "ES", "DEDG", "DEMH", "NO", "DK", "MT", "IT", "CZLU", "PL", "HR", "SK", "LT", "FINW", "ROLS", "EE", "FINM", "FINE", "RODD","TR")) +
  xlab("RDA1 (19.04%)") + ylab("RDA2 (18.08%)") + theme_bw() +
  theme(legend.title = element_text(size = 10),
        legend.position = 'top',
        legend.text = element_text(size = 10)) + theme(panel.grid = element_blank()) + theme(axis.text = element_text(size = 11),
                                                                                      axis.title = element_text(size = 20)) +
  guides(color = guide_legend(nrow = 1, size = 3))
PLOT_LEG

leg <- ggpubr::get_legend(PLOT_LEG)
### Plot the overall 'Figure 2', which combines the RDA, pRDA, and LFMM plots:

library(cowplot)

library(tidyverse)
leg <- plot_grid(leg)
PLOT_AB <- plot_grid(PLOT_A,PLOT_B, labels = 'AUTO', nrow = 1)
PLOT_AB <- plot_grid(leg, PLOT_AB, rel_heights = c(0.10, 0.90), ncol = 1)
IMAGE <- plot_grid(PLOT_AB,PLOT_C, ncol = 1, labels = c('', 'C'))
IMAGE
ggsave("/file/path/to/IMAGE.pdf", plot = IMAGE, device = "pdf", width = 15, height = 10)
