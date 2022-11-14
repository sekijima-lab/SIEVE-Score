#!/bin/sh
#$ -cwd
#$ -l h_node=1
#$ -l h_rt=12:00:00

source ~/.bashrc

PROT=`sed -n ${SGE_TASK_ID}p all.txt`

cp glide-dock_SP_template.in ../${PROT}/glide-dock.in
cd ../${PROT}

sed -i -e "s/__PROTEIN__/${PROT}/g" glide-dock.in

${SCHRODINGER}/glide glide-dock.in -OVERWRITE -NJOBS 10\
    -HOST "localhost:14" -TMPDIR $TMPDIR -ATTACHED -WAIT
