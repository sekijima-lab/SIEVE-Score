#!/bin/sh
#$ -cwd
#$ -l h_node=1
#$ -l h_rt=24:00:00
#$ -o ../logs
#$ -o ../logs

source ~/.bashrc

for PROT in FAK FAK_DUDE; do
    cd ../${PROT}

    ${SCHRODINGER}/run python ../../../SIEVE-Score.py @SIEVE-Score_settings.txt --nprocs 14
done
