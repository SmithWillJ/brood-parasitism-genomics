#!/bin/bash
#SBATCH --job-name=EEMS
#SBATCH --output=EEMS_error.txt
#SBATCH --error=EEMS_out.txt
#SBATCH --account=project_xxxxx
#SBATCH --time=10:00:00

# Load plink2, plus the data

module load plink/2.00a5

plink2 --vcf PLINKPRUNED_Gene_Flow.recode.vcf --allow-extra-chr --make-bed --out PLINKPRUNED_EEMS_input

export PATH=/appl/soft/bio/eems/bin:$PATH

# bed2diffs struggles with non-model chromosomes, so the following reformatting is necessary:

cp PLINKPRUNED_EEMS_input.bed PLINKPRUNED_EEMS_input_new.bed
cp PLINKPRUNED_EEMS_input.fam PLINKPRUNED_EEMS_input_new.fam
awk '{print 1,$2,$3,$4,$5,$6}' PLINKPRUNED_EEMS_input.bim > PLINKPRUNED_EEMS_input_new.bim

bed2diffs_v1 --bfile PLINKPRUNED_EEMS_input_new

#### Export the PLINKPRUNED_EEMS_input_new.diffs file to your Desktop.
#### Edit the settings file (to set number of demes plus output name) with: nano params-simno1.ini
#### Run EEMS with: ./runeems_snps –params params-simno1.ini
#### We will run EEMS eight times, changing the number of demes as follows: 100, 300, 500, 700, 900, 1100, 1300, 1500.
#### We will plot the combined result (i.e. across different variations in the number of demes) in R with the EEMS_Plotting.R script.
