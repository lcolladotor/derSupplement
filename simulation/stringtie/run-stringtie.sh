#!/bin/sh

# Directories
MAINDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/stringtie
DATADIR=${MAINDIR}/hisat

# Define variables
CORES=4

mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

touch ${WDIR}/GTFfiles.txt

ls ${DATADIR}/*[0-9].bam | while read fullpath
    do
    bamfile=$(basename "${fullpath}")
    libname="${bamfile%.*}"
    echo "Creating script for ${libname}"
    sname="${libname}.stringtie"
    
    ## Create GTF file list for cuffmerge
    echo "${WDIR}/${libname}/outfile.gtf" >> ${WDIR}/GTFfiles.txt
    
    ## Create scripts    
	cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=4G,h_fsize=10G
#$ -pe local ${CORES}
#$ -N ${sname}

echo "**** Job starts ****"
date

cd ${WDIR}

## Create output directory
mkdir -p ${WDIR}/${libname}/

## StringTie version used
stringtie --version

## Run StringTie
stringtie ${fullpath} -o ${WDIR}/${libname}/outfile.gtf -p ${CORES} -G ${MAINDIR}/gtf/chr17.gtf

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF
	call="qsub ${WDIR}/.${sname}.sh"
	echo $call
	$call
done


## Now run cuffmerge

sname="cmerge-stie"
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=3G,h_vmem=4G,h_fsize=10G
#$ -pe local ${CORES}
#$ -N ${sname}
#$ -hold_jid *stringtie

echo "**** Job starts ****"
date

cd ${WDIR}

## Load Cuffmerge
module load cufflinks/2.2.1

## Cuffmerge version usedd
cuffmerge --version

## Run cuffmerge
cuffmerge -o ${WDIR}/cuffmerge -p ${CORES} -g ${MAINDIR}/gtf/chr17.gtf ${WDIR}/GTFfiles.txt

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call

ls ${DATADIR}/*[0-9].bam | while read fullpath
    do
    bamfile=$(basename "${fullpath}")
    libname="${bamfile%.*}"
    echo "Creating script for ${libname}"
    sname="${libname}.bgprep"
        
    ## Create scripts    
	cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=4G,h_fsize=10G
#$ -pe local ${CORES}
#$ -N ${sname}
#$ -hold_jid cmerge-stie

echo "**** Job starts ****"
date

cd ${WDIR}

## StringTie version used
stringtie --version

## Run StringTie
stringtie ${fullpath} -o ${WDIR}/${libname}/outfile-run2.gtf -p ${CORES} -G ${WDIR}/cuffmerge/merged.gtf -b ${WDIR}/${libname}/ -e

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF
	call="qsub ${WDIR}/.${sname}.sh"
	echo $call
	$call
done

