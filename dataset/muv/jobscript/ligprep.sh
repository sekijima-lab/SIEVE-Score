#!/bin/sh
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=12:00:00

source ~/.bashrc

PROT=`sed -n ${SGE_TASK_ID}p all.txt`

cd ../${PROT}
$SCHRODINGER/ligprep -ismi ${PROT}.smi -omae ${PROT}_ligands.maegz -WAIT -NJOBS 3 -TMPDIR $TMPDIR
