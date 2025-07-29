setwd("/file/path/to/Desktop/")
 
install.packages('hierfstat')
install.packages('vcfR')
library(hierfstat)
library(vcfR)
 
vcf <- read.vcfR("populations.snps.vcf")
vcf
pop_map <- read.table("LFMM_popmap_no_NA", header=TRUE, stringsAsFactors=TRUE)
pop_map
genind <- vcfR2genind(vcf)
genind
genind@pop <- pop_map$STRATA
genind$pop
hierfstat <- genind2hierfstat(genind)
summary(hierfstat)
out <- write.bayescan(dat=hierfstat,diploid = TRUE, fn='test.bsc')

#### What about the env data (i.e. parasitism rate)?
d <- data.frame(data = c(14.6, 14.5, 8.3, 8.1, 0, 21.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 14.3, 0, 5.8, 0, 0))
d$data <- d$data/sd(d$data)
write.table(d, file = '~/Desktop/test.txt', sep = " ", col.names = F, row.names = F)
