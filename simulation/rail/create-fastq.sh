#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=20G,h_fsize=30G
#$ -N rail-fastq
echo "**** Job starts ****"
date

mkdir -p logs

## Create fake fastq files
module load R/devel
Rscript create-fastq.R

mv rail-fastq.* logs/

echo "**** Job ends ****"
date