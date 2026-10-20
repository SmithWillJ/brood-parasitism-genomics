#!/bin/bash
#### 1. Generate VCF for gene flow and population structure analysis
####    (repeat with different heterozygosity settings)

#SBATCH --job-name=Gene_Flow_VCF_Gen
#SBATCH --output=Gene_Flow_VCF_out.txt
#SBATCH --error=Gene_Flow_VCF_error.txt
#SBATCH --account=project_xxxxx

set -euo pipefail

module load R-uoneasy/4.2.1-foss-2022a      # not used in this script
module load Stacks/2.62-foss-2023a          # ref_map.pl, populations
module load VCFtools/0.1.16-GCC-12.3.0      # VCF stats and filtering
module load BCFtools/1.18-GCC-12.3.0        # genotype masking, sorting, indexing
module load PLINK/1.9b_6.21-x86_64          # LD pruning

cd /filepath/populations_runs

# Output directory for this parameter set

mkdir -p /filepath/populations_runs/refmap_outputs_GEA_auto_0.5

# run the Stacks pipeline

ref_map.pl --samples /filepath/alignments_Will_Trial --popmap /filepath/info/GEA_PopMap -o /filepath/populations_runs/refmap_outputs_GEA_auto_0.5 -X "populations: --vcf"

# One directory per populations run

mkdir -p GEA_auto_0.5_run1   # run1: unfiltered baseline, used for QC stats
mkdir -p GEA_auto_0.5_run2   # run2: depth blacklist applied + HWE stats
mkdir -p GEA_auto_0.5_run3   # run3: final filtered dataset

# populations with no filters - to get a baseline VCF

populations -P /filepath/populations_runs/refmap_outputs_GEA_auto_0.5 --popmap /filepath/info/GEA_PopMap -O /filepath/populations_runs/GEA_auto_0.5_run1 --vcf

cd /filepath/populations_runs/GEA_auto_0.5_run1

# QC statistics

vcftools --vcf populations.snps.vcf --depth --out test-depth                # mean depth per individual
vcftools --vcf populations.snps.vcf --site-mean-depth --out test-mean-depth # mean depth per site
vcftools --vcf populations.snps.vcf --missing-site --out test_missing-site  # missingness per site
vcftools --vcf populations.snps.vcf --missing-indv --out test_missing-indv  # missingness per individual
vcftools --vcf populations.snps.vcf --freq2 --out test_freq                 # allele frequencies

# List sites that FAIL the mean-depth window (mean depth <10 or >100).
# Low depth = unreliable calls; very high depth = likely paralogs/repeats.
# --removed-sites writes the failing sites to file.removed.sites (CHROM, POS).

vcftools --vcf populations.snps.vcf --min-meanDP 10 --max-meanDP 100 --removed-sites --out vcftools_min-meanDP10_max-meanDP100

# Convert the removed sites into a Stacks blacklist of locus IDs:

grep -f vcftools_min-meanDP10_max-meanDP100.removed.sites populations.sumstats.tsv -F -w | cut -f1 | sort -n | uniq > blacklist_min_meanDP10_max-meanDP100_removed_positions

# Re-run populations excluding the depth-blacklisted loci

populations -P /filepath/populations_runs/refmap_outputs_GEA_auto_0.5 --popmap /filepath/info/GEA_PopMap -O /filepath/populations_runs/GEA_auto_0.5_run2 -B blacklist_min_meanDP10_max-meanDP100_removed_positions --vcf --hwe

cd /filepath/populations_runs/GEA_auto_0.5_run2

# build a whitelist of SNPs in HWE (p > 0.05) in the listed populations

cat populations.sumstats.tsv | grep -v '^#' | awk '$20 > 0.05 && ($5 == "DEDG" || $5 == "DEMH" || $5 == "EE" || $5 == "FINM" || $5 == "FINE" || $5 == "IT" || $5 == "LT" || $5 == "NO" || $5 == "PL")' | cut -f 1,4 | awk '++a[$1,$2]==1{ print $1,$2 }' OFS='\t' | sort -n > wl_snps_in_hwe_in_all_pops_with_10indv

# Final populations run with the HWE whitelist and stricter filters

populations -P /filepath/populations_runs/refmap_outputs_GEA_auto_0.5 --popmap /filepath/info/GEA_PopMap -O /filepath/scripts/populations_runs/GEA_auto_0.5_run3/ -W wl_snps_in_hwe_in_all_pops_with_10indv -R 1 --min-mac 3 --max-obs-het 0.5 --write-single-snp --vcf

cd /filepath/populations_runs/GEA_auto_0.5_run3/

# Keep only the listed chromosomes

vcftools --vcf populations.snps.vcf --chr OU383773.1 --chr OU383774.1 --chr OU383775.1 --chr OU383777.1 --chr OU383778.1 --chr OU383779.1 --chr OU383780.1 --chr OU383781.1 --chr OU383782.1 --chr OU383783.1 --chr OU383784.1 --chr OU383785.1 --chr OU383786.1 --chr OU383787.1 --chr OU383788.1 --chr OU383789.1 --chr OU383790.1 --chr OU383791.1 --chr OU383792.1 --chr OU383793.1 --chr OU383794.1 --chr OU383796.1 --chr OU383797.1 --chr OU383798.1 --chr OU383800.1 --chr OU383799.1 --chr OU383801.1 --chr OU383802.1 --chr OU383803.1 --recode --out populations.snps.chromfilt

# ALLELIC BALANCE FILTER
# Sets genotypes to missing (./.) where a binomial test on the allele depths
# (AD) gives p < 1e-3. This removes heterozygous calls with strongly skewed
# read counts (e.g. 18 reads vs 2 reads).
#   -t "b:AD<1e-3" : target genotypes failing the binomial test
#   -n .           : set them to missing

bcftools +setGT -O z -o populations.snps.chromfilt.ABfilt.unsorted.vcf.gz populations.snps.chromfilt.recode.vcf -- -t "b:AD<1e-3" -n .

# Sort and index

bcftools sort populations.snps.chromfilt.ABfilt.unsorted.vcf.gz -O z -o populations.snps.chromfilt.ABfilt.vcf.gz
bcftools index populations.snps.chromfilt.ABfilt.vcf.gz

# Remove sites with >20% missing genotypes (keep sites genotyped in >=80% of samples)

vcftools --gzvcf populations.snps.chromfilt.ABfilt.vcf.gz --max-missing 0.8 --recode --out populations.snps.chromfilt.ABfilt.maxmiss0.8

# LD pruning

plink --vcf populations.snps.chromfilt.ABfilt.maxmiss0.8.recode.vcf --allow-extra-chr --double-id --make-bed --out PLINK_OUTPUT

plink --bfile PLINK_OUTPUT --allow-extra-chr --indep-pairwise 50 5 0.5

# LD-pruned VCF

vcftools --vcf populations.snps.chromfilt.ABfilt.maxmiss0.8.recode.vcf --exclude plink.prune.out --recode --out PLINKPRUNED
