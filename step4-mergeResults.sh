#!/bin/sh

## Usage
# sh step4-mergeResults.sh brainspan run5-v1.5.30

# Define variables
EXPERIMENT=$1
SHORT="derM-${EXPERIMENT}"
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
#$ -l mem_free=80G,h_vmem=100G,h_fsize=10G
#$ -N ${sname}
#$ -hold_jid derA-${EXPERIMENT}.${PREFIX}.chr*

echo "**** Job starts ****"
date

mkdir -p ${WDIR}/${outdir}/logs

# merge results
cd ${WDIR}
module load R/3.3
Rscript -e "library(derfinder); load('/dcl01/lieber/ajaffe/derRuns/derfinderExample/derGenomicState/GenomicState.Hsapiens.UCSC.hg19.knownGene.Rdata'); load('${WDIR}/${PREFIX}/chr22/optionsStats.Rdata'); chrs <- c(1:22, 'X', 'Y'); mergeResults(chrs = chrs, prefix = '${PREFIX}', genomicState = GenomicState.Hsapiens.UCSC.hg19.knownGene[['fullGenome']], optionsStats = optionsStats); Sys.time(); proc.time(); options(width = 120); devtools::session_info()"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
call="qsub .${sname}.sh"
echo $call
$call
