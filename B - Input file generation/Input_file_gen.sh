####1. Generate VCF for gene flow and population structure analysis*

#!/bin/bash
#SBATCH --job-name=Gene_Flow_VCF_Gen
#SBATCH --output=Gene_Flow_VCF_out.txt
#SBATCH --error=Gene_Flow_VCF_error.txt
#SBATCH --account=project_xxxxx

#### First, run ‘ref_map.pl’:

ref_map.pl --samples /file_path_to_list_of_sample_alignments/ --popmap /filepath/to/popmap/ -o /filepath/to/outputs/ -X "populations: --vcf"

cd /filepath/to/populations/runs/

module load biokit
mkdir populations_run1_flow
mkdir populations_run2_flow
mkdir populations_run3_flow

populations -P /filepath/for/refmap/outputs/ --popmap /filepath/to/popmap -O /filepath/to/populations_run1_flow/ --vcf

#### We now have an unfiltered VCF in the ‘populations_run1_flow' folder. Let’s use VCFtools to get some statistics on this VCF:

cd /filepath/to/populations_run1_flow/

module load vcftools

vcftools --vcf populations.snps.vcf --depth --out test-depth
vcftools --vcf populations.snps.vcf --site-mean-depth --out test-mean-depth
vcftools --vcf populations.snps.vcf --missing-site --out test_missing-site
vcftools --vcf populations.snps.vcf --missing-indv --out test_missing-indv
vcftools --vcf populations.snps.vcf --freq2 --out test_freq

#### Now let’s do some depth filtering in VCFtools. First, we should output the sites to be removed:

vcftools --vcf populations.snps.vcf --min-meanDP 10 --max-meanDP 100 --removed-sites --out vcftools_min-meanDP10_max-meanDP100

grep -f vcftools_min-meanDP10_max-meanDP100.removed.sites populations.sumstats.tsv -F -w | cut -f1 | sort -n | uniq > blacklist_min_meanDP10_max-meanDP100_removed_positions

#### Run ‘populations’ with the blacklist we just made, plus some additional filters, and output HWE. We will output the new VCF in ‘populations_run2_flow’:

populations -P /filepath/to/refmap_outputs/ --popmap /filepath/to/popmap -O /filepath/to/populations_run2_flow/ -B /filepath/to/blacklist_min_meanDP10_max-meanDP100_removed_positions --vcf --hwe

#### Next, we want to output a whitelist of SNPs in HWE in all populations which have more than ten individuals. This is done on the new VCF:

cd /filepath/to/populations_run2_flow/

cat populations.sumstats.tsv | grep -v '^#' | awk '$20 > 0.05 && ($5 == "DEDG" || $5 == "DEMH" || $5 == "EE" || $5 == "FINM" || $5 == "FINE" || $5 == "IT" || $5 == "LT" || $5 == "NO" || $5 == "PL")' | cut -f 1,4 | awk '++a[$1,$2]==1{ print $1,$2 }' OFS='\t' | sort -n > wl_snps_in_hwe_in_all_pops_with_10indv

#### Now let’s run populations again with this new whitelist (i.e. filtering our depth filtered VCF). The VCF which is filtered for depth AND hwe AND MAF AND –R AND heterozygosity is outputted to ‘populations_runs3_flow’:

populations -P /filepath/to/refmap_outputs/ --popmap /filepath/to/popmap -O /filepath/to/populations_run3_flow/ -W /filepath/to/wl_snps_in_hwe_in_all_pops_with_10indv -R 1 --min-mac 3 --max-obs-het 0.75 --write-single-snp --vcf

cd /filepath/to/populations_run3_flow/

module load plink/1.90

#### Following this, we use PLINK to filter for linkage disequilibrium – important as this could impact analyses of gene flow or population genetic structure.

module load plink/1.90

plink --vcf populations.snps.vcf --allow-extra-chr --make-bed --out PLINK_OUTPUT
plink --bfile PLINK_OUTPUT --allow-extra-chr --indep-pairwise 50 5 0.5

vcftools --vcf populations.snps.vcf --exclude plink.prune.out --recode --out PLINKPRUNED

#### 2. Generate VCF for gene-environment association analyses*

#!/bin/bash
#SBATCH --job-name=GEA_VCF_Gen
#SBATCH --output=GEA_VCF_out.txt
#SBATCH --error=GEA_VCF_error.txt
#SBATCH --account=project_xxxxx

cd /filepath/to/populations/runs/

module load biokit
mkdir populations_run1_GEA
mkdir populations_run2_GEA
mkdir populations_run3_GEA

populations -P /file/path/torefmap_outputs/ --popmap /file/path/to/lfmm_popmap_no_NA -O /filepath/to/populations_run1_GEA/ --vcf

#### We now have an unfiltered VCF in the ‘populations_run1_GEA’ folder. Let’s use VCFtools to get some statistics on this VCF:

cd /filepath/to/populations_run1_GEA/

module load vcftools

vcftools --vcf populations.snps.vcf --depth --out test-depth
vcftools --vcf populations.snps.vcf --site-mean-depth --out test-mean-depth
vcftools --vcf populations.snps.vcf --missing-site --out test_missing-site
vcftools --vcf populations.snps.vcf --missing-indv --out test_missing-indv
vcftools --vcf populations.snps.vcf --freq2 --out test_freq

#### Now let’s do some depth filtering in VCFtools. First, we should output the sites to be removed:

vcftools --vcf populations.snps.vcf --min-meanDP 10 --max-meanDP 100 --removed-sites --out vcftools_min-meanDP10_max-meanDP100

grep -f vcftools_min-meanDP10_max-meanDP100.removed.sites populations.sumstats.tsv -F -w | cut -f1 | sort -n | uniq > blacklist_min_meanDP10_max-meanDP100_removed_positions

#### Run ‘populations’ with the blacklist we just made, plus some additional filters, and output HWE. We will output the new VCF in ‘populations_run2_GEA_will’:

populations -P /filepath/torefmap_outputs/ --popmap /filepath/to/lfmm_popmap_no_NA -O /filepath/to/populations_run2_GEA/ -B /filepath/to/blacklist_min_meanDP10_max-meanDP100_removed_positions --vcf --hwe

#### Next, we want to output a whitelist of SNPs in HWE in all populations which have more than ten individuals. This is done on the new VCF:

cd /filepath/to/populations_run2_GEA

cat populations.sumstats.tsv | grep -v '^#' | awk '$20 > 0.05 && ($5 == "DEDG" || $5 == "DEMH" || $5 == "EE" || $5 == "FINM" || $5 == "FINE" || $5 == "IT" || $5 == "LT" || $5 == "NO" || $5 == "PL")' | cut -f 1,4 | awk '++a[$1,$2]==1{ print $1,$2 }' OFS='\t' | sort -n > wl_snps_in_hwe_in_all_pops_with_10indv

#### Now let’s run populations again with this new whitelist (i.e. filtering our depth filtered VCF). The VCF which is filtered for depth AND hwe AND MAF AND –R AND heterozygosity is outputted to ‘populations_run3_GEA’:

populations -P /filepath/torefmap_outputs/ --popmap /filepath/to/lfmm_popmap_no_NA -O /filepath/to/populations_run3_GEA/ -W /filepath/to/wl_snps_in_hwe_in_all_pops_with_10indv -R 1 --min-mac 3 --max-obs-het 0.75 --write-single-snp --vcf

cd /filepath/to/populations_run3_GEA_will/

