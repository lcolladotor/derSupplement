#!/bin/bash

#$ -cwd
#$ -l jabba,mem_free=10G,h_vmem=50G,h_fsize=40G 
#$ -N exonsGFF
#$ -m e

echo "**** Job starts ****"
date

## Create GFF files
Rscript createGFF.R

### Done
echo "**** Job ends ****"
date
