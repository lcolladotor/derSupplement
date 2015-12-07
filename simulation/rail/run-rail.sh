#!/bin/sh

# Directories
MAINDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/rail
DATADIR=${MAINDIR}/rail/simulated_fastq

# Define variables
CORES=10
BOWTIE1=/amber2/scratch/jleek/iGenomes-index/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome
BOWTIE2=/amber2/scratch/jleek/iGenomes-index/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs


for replicate in 1 2 3
    do 
    sname="rail-prep-R${replicate}"
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=15G,h_vmem=25G,h_fsize=30G
#$ -pe local ${CORES}
#$ -N ${sname}
#$ -hold_jid rail-fastq
echo "**** Job starts ****"
date

cd ${WDIR}

## run prep
rail-rna --version
rail-rna prep local -m ${WDIR}/rail-manifest-R${replicate}.txt -o sim_prepped_R${replicate} -p ${CORES} --log rail-rna_logs-R${replicate}

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF
    call="qsub ${WDIR}/.${sname}.sh"
    echo $call
    $call
done

for replicate in 1 2 3
    do 
    sname="rail-align-R${replicate}"
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=30G,h_vmem=45G,h_fsize=30G
#$ -pe local ${CORES}
#$ -N ${sname}
#$ -hold_jid rail-prep-R${replicate}
echo "**** Job starts ****"
date

cd ${WDIR}

## run prep
rail-rna --version
rail-rna align local -i sim_prepped_R${replicate} -m ${WDIR}/rail-manifest-R${replicate}.txt -x ${BOWTIE1},${BOWTIE2} -p ${CORES} -o rail-rna_out-R${replicate} --log rail-rna_logs-R${replicate}

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF
    call="qsub ${WDIR}/.${sname}.sh"
    echo $call
    $call
done


