#!/bin/sh

## Usage
# sh step10-region_cuts.sh brainspan

# Define variables
EXPERIMENT=$1
SHORT="regCuts-${EXPERIMENT}"

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derSupplement
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/regionMatrix

if [[ "${EXPERIMENT}" == "brainspan" ]]
then
    echo ""
else
    echo "Specify a valid experiment: brainspan"
fi


# Construct shell files


for chrnum in 22 21 Y 20 19 18 17 16 15 14 13 12 11 10 9 8 X 7 6 5 4 3 2 1
do
    chr="chr${chrnum}"
    sname="${SHORT}-${chr}"
    echo "Creating script ${sname}"

    cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=70G,h_vmem=90G,h_fsize=30G
#$ -N ${sname}
#$ -hold_jid fullCov-${EXPERIMENT}

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load coverage & get regions
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step10-region_cuts.R -m "${MAINDIR}" -c "${chr}"

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
done
    
sname="${SHORT}-merge"

echo "Creating script for merging the region_cuts results"
cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=100G,h_vmem=120G,h_fsize=40G
#$ -N ${sname}
#$ -hold_jid ${SHORT}-chr*

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Combine regions
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step10b-region_cuts-merge.R

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF
    
call="qsub .${sname}.sh"
echo $call
$call
