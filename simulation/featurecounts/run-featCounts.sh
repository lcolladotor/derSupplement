#!/bin/sh

# Directories
MAINDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/featurecounts

# Define variables
CORES=4

mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs


for replicate in 1 2 3
    do 
    sname="featCounts-R${replicate}"
    ## Create script
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=3G,h_vmem=4G,h_fsize=10G
#$ -pe local ${CORES}
#$ -N ${sname}

echo "**** Job starts ****"
date

cd ${WDIR}


## Load R
module load R/3.3

## Run featureCounts
Rscript featureCounts.R -r ${replicate}

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF

    call="qsub ${WDIR}/.${sname}.sh"
    echo $call
    $call
done
