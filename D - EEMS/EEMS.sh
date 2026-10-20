#!/bin/bash
#SBATCH --job-name=eems
#SBATCH --output=eems_%A_%a.out
#SBATCH --error=eems_%A_%a.err
#SBATCH --time=150:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-7

module load R-uoneasy/4.2.1-foss-2022a
module load GCC/10.2.0
module load Boost/1.74.0-GCC-10.2.0
module load Eigen/3.3.8-GCCcore-10.2.0
module load BCFtools/1.18-GCC-12.3.0

VCF=/filepath/populations_runs/GF_auto_0.5_run3/populations.snps.vcf
COORD=/filepath/scripts/Individuals.coords
OUTER=/filepath/scripts/Outer.outer
EEMS=/filepath/eems/runeems_snps/src
OUT=/filepath/eems_runs

mkdir -p $OUT

DEMES=(100 300 500 700 900 1100 1300 1500)
NDEMES=${DEMES[$SLURM_ARRAY_TASK_ID]}
RUN=$OUT/run_${NDEMES}

mkdir -p $RUN
cd $RUN

# extract genotypes

bcftools query -f '[%GT\t]\n' $VCF > geno.txt

# compute diffs

python3 << EOF
import numpy as np

geno = []
with open("geno.txt") as f:
    for line in f:
        g = line.strip().split()
        row = []
        for x in g:
            if x == "./.":
                row.append(np.nan)
            else:
                a,b = x.split("/")
                row.append(int(a)+int(b))
        geno.append(row)

geno = np.array(geno).T

n = geno.shape[0]
D = np.zeros((n,n))

for i in range(n):
    for j in range(i, n):
        valid = ~np.isnan(geno[i]) & ~np.isnan(geno[j])
        if np.sum(valid) == 0:
            d = 0
        else:
            d = np.mean(np.abs(geno[i][valid] - geno[j][valid]))
        D[i,j] = D[j,i] = d

np.savetxt("data.diffs", D, fmt="%.6f")
EOF

# coords + outer

cp $COORD data.coord
cp $OUTER data.outer

# params

cat > params.ini << EOF
datapath = $RUN/data
mcmcpath = $RUN/eems_out

nIndiv = 190
nSites = 17746

nDemes = $NDEMES

numMCMCIter = 10000000
numBurnIter = 20000
numThinIter = 1800
EOF

# run EEMS

$EEMS/runeems_snps --params params.ini
