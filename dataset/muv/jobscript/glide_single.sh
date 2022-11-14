#!/bin/sh
#$ -cwd
#$ -l f_node=1
#$ -l h_rt=12:00:00

source ~/.bashrc

for PROT in FAK FAK_DUDE; do
    
    cd ../${PROT}
    
    ${SCHRODINGER}/glide glide-dock.in -OVERWRITE -NJOBS 48\
        -HOST "localhost:24" -TMPDIR $TMPDIR -ATTACHED -WAIT
done
