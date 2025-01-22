############################################################
############################################################
####### Identify candidate genes from our outlier SNPs #####
######## and then calculate polygenic scores ###############
############################################################
############################################################

library(tidyverse)
setwd("/file/path/to/data")

### Load the VCF.

vcf <- read_table('/file/path/to/populations.snps.vcf', col_names = T) %>%
  dplyr::select(CHROM, POS, ID) %>% mutate(snpID = row_number())

### Load the candidate SNPS (outputs of LFMM, (p)RDA, and BayeScEnv).

candidates <- read_delim('file/path/to/candidate_snps.txt', delim = ' ') %>%
  distinct()
colnames(candidates) <- 'snpID'

### Filter the VCF to keep just the candidate SNPs

snps <- vcf %>% filter(snpID %in% candidates$snpID)

###### Read in the annotation file for the Reed Warbler genome #######
###### The reference genome is available at: https://figshare.com/articles/dataset/Annotation_of_Acrocephalus_scirpaceus_genome/16622302/1 ########

gff <- read_tsv("/file/path/to/acrocephalus_scirpaceus.gff3", skip = 1, col_names = F) %>%
  rename(CHROM = X1)
gff <- left_join(gff,snps, by = 'CHROM')

###### Filter the annotation file to get the genes of interest ######

gff <- gff %>% filter(POS < X5 & POS > X4)
just_genes <- gff %>% filter(X3 == 'gene') %>% separate(X9, into = c('ID', 'Name', 'gene_biotype','gene_id', 'gene_name', 'source_gene', 'source_gene_common_name', 'transcript_modes'), sep = ';')
just_genes <- just_genes %>% filter(str_detect(Name, 'Name=')) %>% dplyr::select(Name, snpID)
just_genes$Name <- gsub('Name=', '', just_genes$Name)
write_delim(just_genes, '~/file/path/just_genes.txt')

##### Get candidate genes positions.
candidate_vcf <- vcf %>% filter(snpID %in% candidates$snpID) %>% write_delim('~/file/path/candidate_vcf.txt')

###################################################
## Now we will calculate polygenic scores #########
###################################################

## 1. Get allele frequencies of outlier loci using PLINK, and get environmental data (THIS IS FOR ALL CANDIDATE LOCI, NOT JUST THOSE WITH A KNOWN GENE)

freq <- read_table('plink.frq.strat') %>%
  filter(CHR %in% candidate_vcf$CHROM & SNP %in% candidate_vcf$ID) %>%
  write_delim('freq.txt')

allele_freq <- read.table("freq.txt", header = T, stringsAsFactors = F)
head(allele_freq)

parasitism_rates <- read.csv("/file/path/to/parasitism_rates.csv", header = T, stringsAsFactors = F)
head(parasitism_rates)

allele_freq$PARASIT <- rep(parasitism_rates$PARASIT, length(unique(allele_freq$SNP)))
allele_freq

allele_freq_parasitism <- allele_freq %>% dplyr::select(SNP, MAF, PARASIT)
head(allele_freq_parasitism)

## 2. perform correlations between allele frequencies at each site, and environmental variation

correlations_by_snp <- allele_freq_parasitism %>%
  group_by(SNP) %>%
  summarize(COR=cor(PARASIT, MAF))

correlations_by_snp

## 3. Convert a vcf with the outlier SNPs to 012 format with vcftools.

## make tab-delimited candidates with chrom and position
tab_cand <- candidate_vcf %>% dplyr::select(CHROM, POS) %>% write_tsv('tab_cand.tsv')

snps <- read.table("/file/path/to/out.012", header = F, stringsAsFactors = F)
dim(snps)
snps <- snps[,-1]

# read in candidate and individual IDs

INDS <- read.table("/file/path/to/out.012.indv", header = F)
INDS

SNPS <- vcf %>% filter(snpID %in% candidates$snpID) %>% dplyr::select(ID)
SNPS

colnames(snps) <- SNPS$ID
rownames(snps) <- INDS$V1

## Invert 0 and 2 for SNPs with a negative correlation 

head(correlations_by_snp)

snps[snps == -1] <- NA
cor_neg <- filter(correlations_by_snp, COR < 0)
cor_neg

snps_neg <- dplyr::select(snps, cor_neg$SNP)
dim(snps_neg)  
for (i in 1:ncol(snps_neg)) {
  snps_neg[,i][which(snps_neg[,i] == 0)] <- 9 
  snps_neg[,i][which(snps_neg[,i] == 2)] <- 0
  snps_neg[,i][which(snps_neg[,i] == 9)] <- 2
}
snps_neg

cor_pos <- filter(correlations_by_snp, COR > 0)
snps_pos <- dplyr::select(snps, cor_pos$SNP)
dim(snps_pos)

snps_recode <- cbind(snps_neg, snps_pos)
dim(snps_recode)
head(snps_recode)

## Calculate polygenic scores

polygenic_score <- as.data.frame(rowSums(snps_recode, na.rm = T))
colnames(polygenic_score) <- "SCORE"
head(polygenic_score)
polygenic_score$IND <- rownames(polygenic_score)

sites <- read_table('/file/path/to/pop.txt', col_names = F) %>% dplyr::select(X1, X2) %>%
  rename(IND = X1, CLST = X2)
polygenic_score <- left_join(polygenic_score, sites, by = 'IND')
head(polygenic_score)

parasitism_sites <- parasitism_rates[, c(2,4)]
parasitism_sites
colnames(parasitism_sites) <- c("CLST", "PARASIT")

polygenic_score <- merge(polygenic_score, parasitism_sites, by = "CLST")
polygenic_score

## Test the relationship between polygenic score and parasitism. 

plot(polygenic_score$PARASIT, polygenic_score$SCORE)

# linear:
mod1 <- lm(polygenic_score$SCORE ~ polygenic_score$PARASIT)
summary(mod1)

# quadratic: 
mod2 <- lm(polygenic_score$SCORE ~ polygenic_score$PARASIT + I(polygenic_score$PARASIT^2))
summary(mod2)

library(MuMIn)

model.sel(mod1,mod2)

model.sel(mod1,mod2,rank = 'BIC')

library(ggplot2)

## Plot  model

(APS_PLOT <- ggplot(data = polygenic_score, aes(x = PARASIT, y = SCORE)) + geom_point() + geom_smooth(method = 'lm') +
    theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size = 20),
                       axis.text = element_text(size = 11)) +
    xlab("Parasitism rate (%)") + ylab("Additive polygenic scores"))

ggsave("/filepath/to/IMAGE.pdf", plot = APS_PLOT_QUAD, device = "pdf", width = 8, height = 5)

### Plot quadratic line

sampleinfo <- read_csv("/file/path/to/sample_info2.csv")

LatitudeInds <- polygenic_score %>% 
  dplyr::select(CLST,IND) %>% 
  rename(Selected = IND) %>% 
  left_join(sampleinfo %>% dplyr::select(Selected, lat), by = 'Selected') %>%
  arrange(lat)

Latcolours <- data.frame(CLST = factor(unique(LatitudeInds$CLST), levels = unique(LatitudeInds$CLST[order(-LatitudeInds$lat)])),
                         colrs = c('deeppink2','red','salmon','orange','brown1','aquamarine','lightgreen',
                                   'brown','green','coral','burlywood','darkorange2','darkgreen','pink','purple','violet', 'darkgrey','blue','darkblue','lightblue'))

polygenic_score_plot <- polygenic_score %>%
  left_join(Latcolours, by = 'CLST')
polygenic_score_plot$CLST <- factor(polygenic_score_plot$CLST, levels = levels(Latcolours$CLST))
color_map <- setNames(polygenic_score_plot$colrs, polygenic_score_plot$CLST)

(APS_PLOT_QUAD <- ggplot(data = polygenic_score_plot, aes(x = PARASIT, y = SCORE)) + geom_jitter(width = 0.6, height = 0.6, aes(col = CLST)) + geom_smooth(method = 'lm', formula = y ~ x + I(x^2)) +
    theme_bw() +  scale_colour_manual(values = color_map, name = 'Location') +
    annotate(x = 0.6, y = 34, label = 'R^2 == 0.68', geom = 'text', parse= T) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 20),
          axis.text = element_text(size = 11)) + xlab("Parasitism rate (%)") + ylab("Additive polygenic scores"))

ggsave('Quadratic.png', path = '~/Desktop/', dpi = 300, width = 8, height = 5, unit = 'in')

############ PLOT THE POLYGENIC SCORES AGAINST OTHER ENVIRONMENTAL VARIABLES (AND REJECTION RATE - rejrate.csv).

sample_info <- read_csv('/file/path/to/sample_info2.csv') %>% rename(IND = Selected)

rejrate <- read_csv('/file/path/to/rejrate.csv', col_select = c('Selected', 'rejrate')) %>% 
  rename(IND = Selected) %>% mutate(across(rejrate, as.numeric))

polygenic_score2 <- polygenic_score %>%
  left_join(sample_info %>% dplyr::select(IND, lat, bio4, bio8, bio13, bio15), by = 'IND') %>%
  left_join(rejrate, by = 'IND') %>%
  pivot_longer(PARASIT:rejrate, names_to = 'var', values_to = 'val')

variables <- c('lat', 'bio4', 'bio8', 'bio13', 'bio15', 'PARASIT', 'rejrate')

#### Test models of each variable in the polygenic score except parasit.
pol_score_2_wide <- polygenic_score2 %>%
  pivot_wider(names_from = 'var', values_from = 'val') # prepare a data.frame in wide format.

### TEST BIO4.
mod3 <- lm(data = pol_score_2_wide, SCORE ~ bio4)
mod4 <- lm (data = pol_score_2_wide, SCORE ~ bio4+ I(bio4^2))
model.sel(mod3,mod4,rank = 'BIC') # indistinguishable, use quadratic
summary(mod4)

test_resid_bio4 <- pol_score_2_wide %>% dplyr::select(bio4) %>% cbind(mod4$residuals)
colnames(test_resid_bio4) <- c('bio4', 'resid')
ggplot(data = test_resid_bio4, aes(x = bio4, y = resid)) + geom_point() + geom_smooth(method = 'lm')

### TEST BIO8
mod5 <- lm(data = pol_score_2_wide, SCORE ~ bio8)
mod6 <- lm (data = pol_score_2_wide, SCORE ~ bio8+ I(bio8^2))
model.sel(mod5,mod6,rank = 'BIC')
summary(mod6)

test_resid_bio8 <- pol_score_2_wide %>% dplyr::select(bio8) %>% cbind(mod6$residuals)
colnames(test_resid_bio8) <- c('bio8', 'resid')
ggplot(data = test_resid_bio8, aes(x = bio8, y = resid)) + geom_point() + geom_smooth(method = 'lm')

### TEST BIO13
mod7 <- lm(data = pol_score_2_wide, SCORE ~ bio13)
mod8 <- lm (data = pol_score_2_wide, SCORE ~ bio13+ I(bio13^2))
model.sel(mod7,mod8,rank = 'BIC')
summary(mod8)

test_resid_bio13 <- pol_score_2_wide %>% dplyr::select(bio13) %>% cbind(mod8$residuals)
colnames(test_resid_bio13) <- c('bio13', 'resid')
ggplot(data = test_resid_bio13, aes(x = bio13, y = resid)) + geom_point() + geom_smooth(method = 'lm')

### TEST BIO15
mod9 <- lm(data = pol_score_2_wide, SCORE ~ bio15)
mod10 <- lm (data = pol_score_2_wide, SCORE ~ bio15+ I(bio15^2))
model.sel(mod9,mod10,rank = 'BIC')
summary(mod10)

test_resid_bio15 <- pol_score_2_wide %>% dplyr::select(bio15) %>% cbind(mod10$residuals)
colnames(test_resid_bio15) <- c('bio15', 'resid')
ggplot(data = test_resid_bio15, aes(x = bio15, y = resid)) + geom_point() + geom_smooth(method = 'lm')

### TEST REJRATE
mod11 <- lm(data = pol_score_2_wide, SCORE ~ rejrate)
mod12 <- lm (data = pol_score_2_wide, SCORE ~ rejrate+ I(rejrate^2))
model.sel(mod11,mod12,rank = 'BIC')
summary(mod12)

test_resid_rejrate <- pol_score_2_wide %>% drop_na(rejrate) %>% dplyr::select(rejrate) %>% cbind(mod12$residuals)
colnames(test_resid_rejrate) <- c('rejrate', 'resid')
ggplot(data = test_resid_rejrate, aes(x = rejrate, y = resid)) + geom_point() + geom_smooth(method = 'lm')

### TEST LATITUDE
mod13 <- lm(data = pol_score_2_wide, SCORE ~ lat)
mod14 <- lm (data = pol_score_2_wide, SCORE ~ lat + I(lat^2))
model.sel(mod13,mod14,rank = 'BIC')
summary(mod13)

test_resid_lat <- pol_score_2_wide %>% drop_na(lat) %>% dplyr::select(lat) %>% cbind(mod13$residuals)
colnames(test_resid_lat) <- c('lat', 'resid')
ggplot(data = test_resid_lat, aes(x = lat, y = resid)) + geom_point() + geom_smooth(method = 'lm')

variables <- c('lat', 'bio4', 'bio8', 'bio13', 'bio15', 'PARASIT', 'rejrate')

(ps1 <- ggplot(data = polygenic_score2 %>% filter(var == 'lat'), 
               aes(x = val, y = SCORE)) + 
    geom_point() + geom_smooth(method = 'lm') +
    theme_bw() + 
    annotate(x = 58, y = 34, label = 'R^2 == 0.040', geom = 'text', parse= T) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 11),
          axis.text = element_text(size = 11)) + xlab('Latitude') + ylab("Add. poly. scores"))

(ps2 <- ggplot(data = polygenic_score2 %>% filter(var == 'bio4'), 
               aes(x = val, y = SCORE)) + 
    geom_point() + geom_smooth(method = 'lm') +
    theme_bw() + 
    annotate(x = 7.5, y = 34, label = 'R^2 == 0.054', geom = 'text', parse= T) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 11),
          axis.text = element_text(size = 11)) + xlab('BIO4') + ylab("Add. poly. scores"))

(ps3 <- ggplot(data = polygenic_score2 %>% filter(var == 'bio8'), 
               aes(x = val, y = SCORE)) + xlim(5, 21) +
    geom_point() + geom_smooth(method = 'lm') +
    theme_bw() + 
    annotate(x = 10, y = 34, label = 'R^2 == 0.050', geom = 'text', parse= T) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 11),
          axis.text = element_text(size = 11)) + xlab('BIO8') + ylab("Add. poly. scores"))

(ps4 <- ggplot(data = polygenic_score2 %>% filter(var == 'bio13'), 
               aes(x = val, y = SCORE)) + 
    geom_point() + geom_smooth(method = 'lm') +
    theme_bw() + 
    annotate(x = 120, y = 34, label = 'R^2 == 0.000', geom = 'text', parse= T) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 11),
          axis.text = element_text(size = 11)) + xlab('BIO13') + ylab("Add. poly. scores"))

(ps5 <- ggplot(data = polygenic_score2 %>% filter(var == 'bio15'), 
               aes(x = val, y = SCORE)) + 
    geom_point() + geom_smooth(method = 'lm') +
    theme_bw() + 
    annotate(x = 75, y = 34, label = 'R^2 == 0.011', geom = 'text', parse= T) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 11),
          axis.text = element_text(size = 11)) + xlab('BIO15') + ylab("Add. poly. scores"))

(ps6 <- ggplot(data = polygenic_score2 %>% filter(var == 'rejrate'), 
               aes(x = val, y = SCORE)) + 
    geom_point() + geom_smooth(method = 'lm') +
    theme_bw() + 
    annotate(x = 0.6, y = 34, label = 'R^2 == 0.034', geom = 'text', parse= T) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 11),
          axis.text = element_text(size = 11)) + xlab('Rejection rate') + ylab("Add. poly. scores"))

plot_grid(ps1, ps2, ps3, ps4, ps5, ps6, ncol = 2, 
          labels = c('A)', 'B)', 'C)', 'D)', 'E)', 'F)'), scale = 0.9)
ggsave('Composite_plot.png', path = '~/output/path/', dpi = 300, width = 7, height = 8, unit = 'in')

##### Where is Ase64 microsatellite? We do not have a SNP at this site.

## From Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi), Ase64 is at positions 22090624 - 22089985
mid_micro = (22090624 + 22089985)/2

chrom_8_snps <- vcf %>% filter(CHROM == 'SUPER_8') %>%
  mutate(distance_to_micro = abs(POS - mid_micro),
         distance_to_start = abs(POS - 22089985),
         distance_to_end = abs(POS - 22090624))

ggplot(data = chrom_8_snps, aes(x = distance_to_micro)) + geom_histogram(bins = 100)
