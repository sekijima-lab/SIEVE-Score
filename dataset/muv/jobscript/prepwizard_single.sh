#!/bin/sh
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=12:00:00

source ~/.bashrc 

PROT=FAK
PDB=4Q9S
cd ../${PROT}

$SCHRODINGER/utilities/prepwizard -WAIT\
        -disulfides -fillsidechains -fillloops\
        ${PDB}.pdb ${PROT}_protein.maegz
