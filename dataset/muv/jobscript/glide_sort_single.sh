#!/bin/sh
#$ -cwd
#$ -l s_core=1
#$ -l h_rt=6:00:00
#$ -o ../logs
#$ -e ../logs

source ~/.bashrc

for PROT in FAK FAK_DUDE; do
    cd ../${PROT}
    
    ${SCHRODINGER}/utilities/glide_sort -best_by_title \
	-o glide-dock_best_pv.maegz -r glide-dock_best.rept glide-dock_pv.maegz
done
