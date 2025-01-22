#!/bin/bash
#SBATCH --job-name=ADMIXTURE
#SBATCH --output=ADMIXTURE_error.txt
#SBATCH --error=ADMIXTURE_out.txt
#SBATCH --account=project_xxxxx

module load admixture

for i in {1..22}
do
 admixture --cv PLINKPRUNED_EEMS_input_new.bed $i > log${i}.out
done

