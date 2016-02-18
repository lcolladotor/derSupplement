#!/bin/sh

## Usage
# sh optional2-HTSeq.sh snyder
# sh optional2-HTSeq.sh hippo

# Define variables
EXPERIMENT=$1
SHORT="HTSeq-${EXPERIMENT}"

# Directories
ROOTDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/HTSeq

if [[ "${EXPERIMENT}" == "snyder" ]]
then
    SAMFILES='/dcs01/ajaffe/Snyder/RNAseq/TopHat/*out'
elif [[ "${EXPERIMENT}" == "hippo" ]]
then
    SAMFILES='/dcs01/ajaffe/Hippo/TopHat/*out'
else
    echo "Specify a valid experiment: snyder or hippo"
fi


# Construct shell files
for sam in $SAMFILES
do
	current=${sam##*/}
    sname="${SHORT}.${current}"
	echo "Creating script ${sname}"
	cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=50G,h_vmem=100G,h_fsize=30G
#$ -hold_jid sortSam-${EXPERIMENT}
#$ -N ${sname}

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Copy files to a scratch disk
TDIR=\${TMPDIR}/htseq
mkdir -p \${TDIR}
cd \${TDIR}

# Specify file names
orifile=${sam}/accepted_hits_sorted.sam
newfile=\${TDIR}/accepted_hits_sorted.sam
	
# Actually do the copying
echo "Copying \${orifile} to scratch disk \${newfile}"
cp \${orifile} \${newfile}


# run htseq
echo "**** HTSeq starts ****"
date

htseq-count --stranded=no --type='sequence_feature' --idattr='ID' accepted_hits_sorted.sam ${ROOTDIR}/gff/GenomicState.Hsapiens.ensembl.GRCh37.p11.exons.gff > htseq_output.txt

echo "**** HTSeq finishes ****"
date

# Copy back the results
mv htseq_output.txt ${sam}/

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

    call="qsub .${sname}.sh"
    echo $call
    $call
done
