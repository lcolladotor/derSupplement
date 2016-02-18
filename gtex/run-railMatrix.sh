#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=4G,h_vmem=5G,h_fsize=5G
#$ -N rmat-gtex36
#$ -pe local 8

## You need to run create_meanCov.R first!

echo '**** Job starts ****'
date

# Make logs directory
mkdir -p /dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/logs

# Load the data, save the coverage without filtering, then save each file separately
cd /dcs01/ajaffe/Brain/derRuns/derSupplement/gtex
module load R/3.2.x
Rscript railMatrix.R

## Move log files into the logs directory
mv /dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/rmat-gtex36.* /dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/logs/

echo '**** Job ends ****'
date
