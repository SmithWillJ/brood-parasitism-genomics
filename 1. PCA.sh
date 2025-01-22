#!/bin/bash
#SBATCH --job-name=PCA
#SBATCH --output=PCA_error.txt
#SBATCH --error=PCA_out.txt
#SBATCH --account=project_xxxx
#SBATCH --time=05:00:00

plink --vcf PLINKPRUNED_Gene_Flow.recode.vcf --allow-extra-chr --make-bed --out PLINKPRUNED_Gene_Flow_pca_input

plink --bfile PLINKPRUNED_Gene_Flow_pca_input --allow-extra-chr --pca --out PCA_PLINKPRUNED
