#!/bin/sh

# Directories
MAINDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/stringtie
DATADIR=${MAINDIR}/hisat

# Define variables
CORES=4

# GTF file downloaded from
# http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=456818617_6W2Sal1hUnwDVyyPrV9tyX3djUy4&clade=mammal&org=&db=hg19&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=

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
#$ -l mem_free=3G,h_vmem=15G,h_fsize=30G
#$ -pe local ${CORES}
#$ -N ${sname}

echo "**** Job starts ****"
date

cd ${WDIR}

## StringTie version used
stringtie --version

## Run StringTie
stringtie ${fullpath} -o ${WDIR}/${libname}/outfile.gtf -p ${CORES} -G ${WDIR}/hg19-knownGene.GTF

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
#$ -l mem_free=3G,h_vmem=15G,h_fsize=30G
#$ -pe local ${CORES}
#$ -N ${sname}
#$ -hold_jid *stringtie

echo "**** Job starts ****"
date

cd ${WDIR}

## Load Cuffmerge
module avail cufflinks/2.2.1

## Cuffmerge version usedd
cuffmerge --version

## Run cuffmerge
cuffmerge -o ${WDIR}/cmerge-summary -p ${CORES} -g ${WDIR}/hg19-knownGene.GTF ${WDIR}/GTFfiles.txt

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF


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
#$ -l mem_free=3G,h_vmem=15G,h_fsize=30G
#$ -pe local ${CORES}
#$ -N ${sname}
#$ -hold_jid cmerge-stie

echo "**** Job starts ****"
date

cd ${WDIR}

## StringTie version used
stringtie --version

## Run StringTie
stringtie ${fullpath} -o ${WDIR}/${libname}/outfile-run2.gtf -p ${CORES} -G ${WDIR}/merged.gtf -b ${WDIR}/${libname}/ -e

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF
	call="qsub ${WDIR}/.${sname}.sh"
	echo $call
	$call
done

