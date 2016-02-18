#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=4G,h_vmem=5G,h_fsize=5G
#$ -N rmat-gtex36
#$ -pe local 10

echo '**** Job starts ****'
date

# Make logs directory
mkdir -p /dcs01/ajaffe/Brain/derRuns/railDER/gtex36/logs

# Load the data, save the coverage without filtering, then save each file separately
cd /dcs01/ajaffe/Brain/derRuns/railDER/gtex36
module load R/3.2.x
Rscript railMatrix.R

## Move log files into the logs directory
mv /dcs01/ajaffe/Brain/derRuns/railDER/gtex36/rmat-gtex36.* /dcs01/ajaffe/Brain/derRuns/railDER/gtex36/logs/

echo '**** Job ends ****'
date
