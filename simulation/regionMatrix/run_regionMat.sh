#!/bin/sh

# Directories
MAINDIR=/dcl01/lieber/ajaffe/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/regionMatrix

# Define variables
CORES=4

# Create logs dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

for replicate in 1 2 3
    do
    sname="make-regMat-R${replicate}"
    ## Create script
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=3G,h_vmem=4G,h_fsize=10G
#$ -N ${sname}
#$ -pe local ${CORES}

echo "**** Job starts ****"
date

cd ${WDIR}

## Create fullCov
Rscript regionMat.R -r ${replicate}

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF

    call="qsub ${WDIR}/.${sname}.sh"
    echo $call
    $call
done
    