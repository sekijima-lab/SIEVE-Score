#!/bin/bash

qsub -g tga-science -t 1-7 -tc 7 ./glide_sort.sh
