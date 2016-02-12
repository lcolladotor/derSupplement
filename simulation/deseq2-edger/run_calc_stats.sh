#!/bin/sh

# Directories
MAINDIR=/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/deseq2-edger

# Create logs dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

for replicate in 1 2 3
    do
    for complete in yes no
        do 
        sname="stats-featCount-R${replicate}-${complete}"
        ## Create script
        cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=3G,h_vmem=4G,h_fsize=10G
#$ -N ${sname}
#$ -hold_jid featCounts-R${replicate}

echo "**** Job starts ****"
date

cd ${WDIR}

## Calculate stats
Rscript calc_stats.R -r ${replicate} -c "${complete}" -p "featureCounts"

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF

        call="qsub ${WDIR}/.${sname}.sh"
        echo $call
        $call
    done
    
    for pipeline in regionMatrix railMatrix
        do 
        sname="stats-${pipeline}-R${replicate}"
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

## Calculate stats
Rscript calc_stats.R -r ${replicate} -p "${pipeline}"

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF

        call="qsub ${WDIR}/.${sname}.sh"
        echo $call
        $call
    done
done
    