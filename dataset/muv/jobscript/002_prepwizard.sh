#!/bin/sh
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=12:00:00

source ~/.bashrc 

for i in `seq 1 7`; do

    PROT=`sed -n ${i}p all.txt`
    PDB=`sed -n ${i}p pdb.txt`
    cd ../${PROT}
    
    wget https://files.rcsb.org/download/${PDB}.pdb
    
    $SCHRODINGER/utilities/prepwizard\
        -disulfides -fillsidechains -fillloops\
        ${PDB}.pdb ${PROT}_protein.maegz
done
