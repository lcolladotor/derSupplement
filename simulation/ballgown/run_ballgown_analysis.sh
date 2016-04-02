#!/bin/sh

# Directories
MAINDIR=/dcl01/lieber/ajaffe/derRuns/derSupplement/simulation
WDIR=${MAINDIR}/ballgown

# Create logs dir
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

## Run ballgown
Rscript ballgown_analysis.R -r ${replicate} -c "${complete}"

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date

EOF

        call="qsub ${WDIR}/.${sname}.sh"
        echo $call
        $call
    done
done
    