#!/bin/sh
#$ -cwd
#$ -l h_node=1
#$ -l h_rt=24:00:00
#$ -o ../logs
#$ -o ../logs

source ~/.bashrc

PROT=`sed -n ${SGE_TASK_ID}p all.txt`
cd ../${PROT}

${SCHRODINGER}/run python ../../../SIEVE-Score.py @SIEVE-Score_settings.txt --nprocs 14
