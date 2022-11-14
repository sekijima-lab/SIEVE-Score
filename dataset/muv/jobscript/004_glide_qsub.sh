#!/bin/bash

qsub -g tga-science -t 1-7 -tc 1 ./glide.sh
