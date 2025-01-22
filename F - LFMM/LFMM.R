### Parasitism_analyses
### 2024 - WJ Smith and K Rönkä
### Gene-environment association analysis - LFMM
### We will use climate variables and parasitism risk to look for signs of local adaptation.
### Load required libraries and functions
.libPaths(c("/projappl/project_xxxxx/Rpackages_lfmm_4.2.1", .libPaths()))
#BiocManager::install("LEA")

require(raster)
require(gdm)
require(tidyverse)
require(lfmm)
require(adegenet)
require(hierfstat)
require(qqman)
require(rgdal)
require(SimPhe)
require(LEA)

# function to get candidates
get_candidates <- function(pval, fdr){
  pval <- pval[!is.na(pval)]
  L <- length(pval)
  #Benjamini-Hochberg algorithm for FDR
  w <- which(sort(pval) < fdr * (1:L) / L)
  candidates <- order(pval)[w]
  candidates
}

#function to remove intercept from gdm predictions
removeIntercept <- function(mod,pred){
  adjust <- 0 - log(1-pred) - mod$intercept
  adjustDissim <- 1-exp(0-adjust)
  return(adjustDissim)
}

# function to convert allele counts to snp genotype
counts2geno <- function(x){
  ifelse(x==0,"0:0",
         ifelse(x==1, "0:1","1:1"))
}
#### Import and prepare data
### genetic data
#### First let’s import the genomic data and check it out.
snp_matrix <- read.geno("/file/path/to/populations.snps.geno")

dim(snp_matrix)
### Each row is an individual tree, and each column contains allele counts for a single SNP locus. The allele counts are in the format 0/1/2 if the individual has 0, 1, or 2 copies of the ancestral allele.
#### population and climate data
### Now we read in a dataframe giving the population and geographic coordinates of each sample
sample_info2 <- read.csv("/filepath/to/RW_LFMM_sampleinfo2_prec.csv", header = TRUE, sep = ";", dec = ".")
sample_info2$parasit <- as.numeric(gsub(",",".",sample_info2$parasit))
head(sample_info2)
##### For each population, we will extract bioclimatic variables that were downloaded from worldclim.org. First we define the variables we want, and then read them in. The data are in a geospatial raster format. First take in all variables, get values and then test for correlations.
predNames <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
presClim <- stack(paste("/scratch/project_2001434/2021+2017RAD_91bp/analyses/LEA_analyses/wc2.1_30s_bio/", predNames, ".tif", sep=""))
#### Plot the bio4 layer and our population locations to see what it looks like. Try plotting the variables…
plot(presClim[["bio1"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio2"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio3"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio4"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio5"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio6"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio7"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio8"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")], pch=sample_info2$parasit_yes)

plot(presClim[["bio9"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio10"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio11"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio12"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio13"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio14"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio15"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio16"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio17"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio18"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])

plot(presClim[["bio19"]], xlim=c(-15,45), ylim=c(35,68))
points(sample_info2[,c("lon", "lat")])
 #### Now we will use the geographic coordinates in our sample_info dataframe to extract values of the bioclimatic variables for each population. This might take a minute to run
#autosomes
sample_info2 <- cbind(sample_info2,raster::extract(presClim, y=sample_info2[,c("lon","lat")])) %>%
  mutate(bio4=bio4/100) # bio4 is units degrees Celcius x 100, so here we divide by 100 to get back into C 

write.csv(sample_info2, 'sample_info2.csv')

#testing for correlations between the variables
require(GGally) #default is pairwise Pearson correlations
ggcorr(sample_info2[,3:24]) 

cor(sample_info2[,3:24])
ggcorr(sample_info2[c("parasit","lat","bio4","bio8","bio13","bio15")])

cor(sample_info2[c("parasit","lat","bio4","bio8","bio13","bio15")])
alldata <- cbind(sample_info2, snp_matrix)
head(alldata)
sum(is.na(alldata))
sample_info2 <- alldata[,1:24]
predNames <- c("parasit","lat","bio4","bio8","bio13","bio15")
sample_info2[, predNames] <- scale(sample_info2[, predNames])
head(sample_info2)
#testing for correlations between the variables
require(GGally) #default is pairwise Pearson correlations
ggcorr(sample_info2[,3:24]) 

ggcorr(sample_info2[c("parasit","lat","bio4","bio8","bio13","bio15")])

cor(sample_info2[c("parasit","lat","bio4","bio8","bio13","bio15")])
#### Run LFMM to identify climate-related SNPs
##### Now that we have imported and prepped our data, we are ready to run lfmm2 to identify SNPs strongly associated with our climate variables. This program runs a regression model for each SNP. In each model, the SNP is the independent response variable, the climate variables are the dependent predictors, and latent factors is included to control for confounding population structure. Population structure describes allele frequency differences among populations that occur due to neutral processes.
#autosomes
pc <- prcomp(snp_matrix, scale. = TRUE)
screeplot(pc)
##### In our data, we should use K=1. Now we run the lfmm2. The first part (lfmm_ridge) estimates the model, and the second part (lfmm_test) calculates p-values for each of the SNPs. We print out the genomic inflation factors (GIF) for each model. The GIF should be around 1 if population structure is adequately controlled for. If GIF is less than 1, the model is conservative and might missidentify true candidate SNPs as non-significant. If GIF is greater than 1 the model is liberal and we might see many false positives in our results (i.e. SNPs identified as signficant candidates but in reality are not associated with climate adaptation)
#autosomes, using lfmm2
mod.lfmm <- lfmm_ridge(Y = snp_matrix, X = sample_info2[,predNames], K = 1)
mod_test <- lfmm_test(Y = snp_matrix, X = sample_info2[,predNames], lfmm=mod.lfmm)
mod_test$gif
### In the next step, we adjust the p-values by the GIF values to make our p-value histograms more uniform. This is done automatically by lfmm2 and found under $calibrated.pvalue of the model result. We can check that the p-values are well calibrated after this adjustment by looking at their distributions in qq-plots. We should see that most values follow the uniform-distribution line, with a few outliers that contain our candidates.
adjusted_p <- mod_test$calibrated.pvalue
head(adjusted_p)
str(adjusted_p)
par(mfrow=(c(2,4)))
qq(adjusted_p[,"parasit"])
qq(adjusted_p[,"lat"])
qq(adjusted_p[,"bio4"])
qq(adjusted_p[,"bio8"])
qq(adjusted_p[,"bio13"])
qq(adjusted_p[,"bio15"])
#### We adjust for multiple testing using the Benjamin-Hochberg Procedure and tolerate a false discovery rate of 1%. The function below returns candidate loci that remain significant after this adjustment
candidates <- apply(adjusted_p, 2, get_candidates, fdr=0.1)
print(candidates)
#### How many candidates were found for each variable?  Now we visualise our candidates with manhattan plots. On the x-axis is the position of each SNP in each chromosome along the genome, and on the y-axis is the -log(10) of p-values. Our candidates, are highlighted in green.
# Check the number of rows in the SNP position file
str(adjusted_p)
summary(adjusted_p)
snp_positions <- read.table("/scratch/project_2001434/retrieved_2021+2017RAD_91bp/2021+2017RAD_91bp/populations_runs/populations_run3_GEA_will/POS_FIL_EDITED.012.pos")
str(snp_positions)
manhattan_df <- cbind(snp_positions,adjusted_p)
str(manhattan_df)
manhattan_df$V4 <- rownames(manhattan_df)

library(tidyverse)
write_csv(manhattan_df, "ooglie_booglies.csv")

# Manhattan plot
model <- "parasit"
model
candidates[[model]]
par(mfrow=c(1,1))
manhattan(manhattan_df,
          chr="V3", bp="V2", snp="V4", p=model,
          suggestiveline = T, genomewideline = T,
          highlight = candidates[[model]])
