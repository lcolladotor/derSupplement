#!/bin/sh

## Usage
# sh step8-coverageToExon.sh brainspan

# Define variables
EXPERIMENT=$1
SHORT="covToEx-${EXPERIMENT}"
CORES=1

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/Brain/derRuns/derSupplement
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/coverageToExon

if [[ "${EXPERIMENT}" == "brainspan" ]]
then
    RLENGTH=100
else
    echo "Specify a valid experiment: brainspan"
fi


for anno in ensembl ucsc
do
    # Construct shell files
    sname="${SHORT}-${anno}"
    echo "Creating script ${sname}"

    cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=240G,h_vmem=270G,h_fsize=30G
#$ -N ${sname}
#$ -hold_jid fullCov-${EXPERIMENT}

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Run coverageToExon()
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step8-coverageToExon.R -e "${EXPERIMENT}" -a "${anno}" -r ${RLENGTH} -c ${CORES}

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
done
