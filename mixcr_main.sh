#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mem_free=12G,hostname=hpc03

input_dir=fastq_final/
mixcr_results=fastq_final_results/

sample_script=scripts/sample_fasta.R
sample_info=sample_info_final.txt

mkdir $mixcr_results

./mixcr1.sh $input_dir $mixcr_results
seed=1
#Rscript $sample_script $mixcr_results IGH $sample_info $seed
#Rscript $sample_script $mixcr_results IGK $sample_info $seed
#Rscript $sample_script $mixcr_results IGL $sample_info $seed
./mixcr_seed.sh $seed $mixcr_results $mixcr_results
