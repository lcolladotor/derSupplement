## Usage
# sh step7-regMatVsDERs.sh brainspan run5-v1.5.30

# Define variables
EXPERIMENT=$1
PREFIX=$2
SHORT="regVsDERs-${EXPERIMENT}"
ncore=5
cores="${ncore}cores"

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derSupplement
MAINDIR=${ROOTDIR}/${EXPERIMENT}


if [[ "${EXPERIMENT}" == "brainspan" ]]
then
    CUTOFFS="0.1 0.25"
else
    echo "Specify a valid experiment: brainspan"
fi

# Construct shell files
for CUTOFF in ${CUTOFFS}
do
    sname="${SHORT}.${PREFIX}-cut-${CUTOFF}"
    echo "Creating script ${sname}"

    WDIR=${MAINDIR}/regionMatrix-vs-DERs/cut${CUTOFF}-vs-${PREFIX}

cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=200G,h_vmem=300G,h_fsize=30G
#$ -N ${sname}
#$ -hold_jid regMat-${EXPERIMENT}-cut-${CUTOFF}-merge,derM-${EXPERIMENT}.${PREFIX}
echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Compare DERs vs regionMatrix
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step7-regMatVsDERs.R -m "${MAINDIR}" -r "${PREFIX}" -w "${WDIR}" -d "${ROOTDIR}" -c ${CUTOFF}

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
done
