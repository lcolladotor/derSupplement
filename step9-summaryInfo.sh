## Usage
# sh step9-summaryInfo.sh brainspan run5-v1.5.30

# Define variables
EXPERIMENT=$1
PREFIX=$2
SHORT="summInfo-${EXPERIMENT}"

# Directories
ROOTDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement
MAINDIR=${ROOTDIR}/${EXPERIMENT}

# Construct shell files
sname="${SHORT}.${PREFIX}"
echo "Creating script ${sname}"

if [[ "${EXPERIMENT}" == "brainspan" ]]
then
    EXAMPLES='c("the complexity induced by alternative transcription" = 5, "coverage dips" = 16, "and coverage variability even on long single exon regions" = 18)'
else
    echo "Specify a valid experiment: brainspan"
fi

WDIR=${MAINDIR}/summaryInfo/${PREFIX}

cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=130G,h_vmem=160G,h_fsize=10G
#$ -N ${sname}
#$ -hold_jid derM-${EXPERIMENT}.${PREFIX}
echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Compare DERs vs regionMatrix
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step9-summaryInfo.R -s '${EXPERIMENT}' -r '${PREFIX}' -p '${EXAMPLES}' -v TRUE

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

call="qsub .${sname}.sh"
echo $call
$call
