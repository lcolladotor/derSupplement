#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=30G,h_vmem=100G
#$ -N PNAS-hippo
echo "**** Job starts ****"
date

# Generate HTML
Rscript -e "library(rmarkdown); render('compareVsPNAS.Rmd', clean = FALSE)"

echo "**** Job ends ****"
date
