#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=20G,h_vmem=30G,h_fsize=20G
#$ -N render-genReads
echo "**** Job starts ****"
date

# Generate HTML
module load R/devel
Rscript -e "rmarkdown::render('generateReads.Rmd')"


echo "**** Job ends ****"
date
