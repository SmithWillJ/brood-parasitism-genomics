#!/bin/bash
#SBATCH --job-name=process_radtags_script
#SBATCH --output=process_radtags_script_out.txt
#SBATCH --error=process_radtags_script_error.txt
#SBATCH --account=project_xxxx

module load biokit

cd ./

process_radtags -f /filepath/to/sequence_archive_lane1.fastq.gz -o /filepath/to/2017_lane1/ -b /filepath/to/2017_barcodes -e sbfI -r -c -q -t 91

process_radtags -f /filepath/to/sequence_archive_lane2.fastq.gz -o /filepath/to/2017_lane2/ -b /filepath/to/2017_barcodes -e sbfI -r -c -q -t 91

process_radtags -f /filepath/to/2021RADsequences_plate1_lane1.fastq.gz -o /filepath/to/2021_plate1_lane1/ -b /filepath/to/2021_plate1_barcodes -e sbfI -r -c -q -t 91

process_radtags -f /filepath/to/2021RADsequences_plate1_lane2.fastq.gz -o /filepath/to/2021_plate1_lane2/ -b /filepath/to/2021_plate1_barcodes -e sbfI -r -c -q -t 91

process_radtags -f /filepath/to/2021RADsequences_plate2_lane1.fastq.gz -o /filepath/to/2021_plate2_lane1/ -b /filepath/to/2021_plate2_barcodes -e sbfI -r -c -q -t 91

process_radtags -f /filepath/to/2021RADsequences_plate2_lane2.fastq.gz -o /filepath/to/2021_plate2_lane2/ -b /filepath/to/2021_plate2_barcodes -e sbfI -r -c -q -t 91

### The above script rescues barcodes and RAD-tag cutsites (-r), cleans data by removing any read with an uncalled base (-c), discards reads with a low quality score (-q), and truncates the final read length to 91 (-t 91).

### We need to create a list of sample names for each plate. To do so, go into each of the folders of raw data, and run the relevant command:

ls | sed -n 's/\.fq.gz$//p' > namelist_2017
ls | sed -n 's/\.fq.gz$//p' > namelist_2021_plate1
ls | sed -n 's/\.fq.gz$//p' > namelist_2021_plate2

The three output files live here: 

Combine:

for s in `cat /filepath/to/namelist_2021_plate1` ; do cat /filepath/to/2021_plate1_lane1/${s}.fq.gz /filepath/to/2021_plate1_lane2/${s}.fq.gz >> /filepath/to/2021_plate1_combined/$s.fq.gz ; done

for s in `cat /filepath/to/namelist_2021_plate2` ; do cat /filepath/to/2021_plate2_lane1/${s}.fq.gz /filepath/to/2021_plate2_lane2/${s}.fq.gz >> /filepath/to/2021_plate2_combined/$s.fq.gz ; done

for s in `cat /filepath/to/namelist_2017` ; do cat /filepath/to/2017_lane1/${s}.fq.gz /filepath/to/2017_lane2/${s}.fq.gz >> /filepath/to/2017_combined/$s.fq.gz ; done

#### Next, we need to unzip and index the reference genome (see the relevant paper for this reference genome at: Lo Cascio Sætre et al. 2021, Genome Biology and Evolution):

gunzip -c bAcrSci1_1.20210512.curated_primary.fa.gz > bAcrSci1_1.20210512.curated_primary.fasta
bwa index ../genome/bAcrSci1_1.20210512.curated_primary.fasta

### Create a list of all sample names:

ls | sed -n 's/\.fq.gz$//p' > namelist_2021+2017

#### Align to the reference genome (this requires all the samples to be in one folder):

module load bwa/0.7.17

name=$(sed -n ${SLURM_ARRAY_TASK_ID}p /filepath/to/namelist_2021+2017)

bwa mem -M /filepath/to/bAcrSci1_1.20210512.curated_primary.fasta \
/filepath/to/all_samples_combined/${name}.fq.gz | samtools sort -o /filepath/to/outputs/${name}.bam

### We now have per-individual ‘bam’ alignment files.
