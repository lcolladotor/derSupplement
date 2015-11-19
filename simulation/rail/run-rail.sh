#!/bin/sh

# Directories
MAINDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/rail
DATADIR=${MAINDIR}/simulated_reads

# Define variables
CORES=4
BOWTIE1=/amber2/scratch/jleek/iGenomes-index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome
BOWTIE2=/amber2/scratch/jleek/iGenomes-index/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs


cat > ${WDIR}/.rail-prep.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=3G,h_vmem=15G,h_fsize=30G
#$ -N rail-prep
echo "**** Job starts ****"
date

cd ${WDIR}

## run prep
rail-rna --version
rail-rna prep local -m ${WDIR}/rail-manifest.txt -o sim_prepped -p 1

mv ${WDIR}/rail-prep.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF
call="qsub ${WDIR}/.rail-prep.sh"
echo $call
$call



cat > ${WDIR}/.rail-align.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=3G,h_vmem=15G,h_fsize=30G
#$ -pe local ${CORES}
#$ -N rail-align
#$ -hold_jid rail-prep
echo "**** Job starts ****"
date

cd ${WDIR}

## run prep
rail-rna --version
rail-rna align local -i sim_prepped -m ${WDIR}/rail-manifest.txt -x ${BOWTIE1},${BOWTIE2} -p ${CORES}

mv ${WDIR}/rail-align.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF
call="qsub ${WDIR}/.rail-align.sh"
echo $call
$call


