## Usage
# sh step7-regMatVsDERs.sh brainspan run5-v1.5.30

# Define variables
EXPERIMENT=$1
PREFIX=$2
SHORT="regVsDERs-${EXPERIMENT}"
ncore=5
cores="${ncore}cores"

# Directories
ROOTDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement
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
#$ -hold_jid regMat-${EXPERIMENT},derM-${EXPERIMENT}.${PREFIX}
echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Compare DERs vs regionMatrix
cd ${WDIR}
module load R/devel
Rscript -e "analysisPath <- '${WDIR}'; load('${MAINDIR}/regionMatrix/regionMat-cut${CUTOFF}.Rdata'); proc.time(); load('${MAINDIR}/derAnalysis/${PREFIX}/fullRegions.Rdata'); proc.time(); library(rmarkdown); render('${ROOTDIR}/step7-regMatVsDERs.Rmd', output_file='${WDIR}/step7-regMatVsDERs.html')"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
done
