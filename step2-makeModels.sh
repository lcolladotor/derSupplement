#!/bin/sh

## Usage
# sh step2-makeModels.sh brainspan run5-v1.5.30

# Define variables
EXPERIMENT=$1
SHORT="derMod-${EXPERIMENT}"
PREFIX=$2

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derSupplement
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/derAnalysis

# Construct shell files
outdir="${PREFIX}"
sname="${SHORT}.${PREFIX}"
echo "Creating script ${sname}"
cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=5G,h_vmem=7G,h_fsize=10G
#$ -N ${sname}
#$ -hold_jid fullCov-${EXPERIMENT}

echo "**** Job starts ****"
date

mkdir -p ${WDIR}/${outdir}/logs

# merge results
cd ${WDIR}/${outdir}/
module load R/3.3
Rscript ${ROOTDIR}/step2-makeModels.R -e "${EXPERIMENT}"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
call="qsub .${sname}.sh"
echo $call
$call
