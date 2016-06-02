#!/bin/bash

filepath=~/ACLS/DUD-E

for i in `cat div.txt`; do
    mkdir $i
    cd $i

    cp $filepath/$i/docking/DUD_dock_merged_best_pv.maegz ./docking_pv.maegz
    cp -r $filepath/$i/SIFt .
    cp $filepath/$i/docking/actives.txt .
    cp $filepath/$i/*_final.ism .
    cp $filepath/$i/*_final.sdf .
    cd ..
    echo "$i end."
done