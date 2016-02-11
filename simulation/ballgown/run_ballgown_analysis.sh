#!/bin/sh

# Directories
MAINDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/ballgown

# Create log dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

for replicate in 1 2 3
    do
    for complete in yes no
        do 
        sname="bg-R${replicate}-${complete}"
        ## Create script
        cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=3G,h_vmem=4G,h_fsize=10G
#$ -N ${sname}

echo "**** Job starts ****"
date

cd ${WDIR}

## Run featureCounts
Rscript ballgown-analysis.R -r ${replicate} -c "${complete}"

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF

        call="qsub ${WDIR}/.${sname}.sh"
        echo $call
        $call
    done
done
    