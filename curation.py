#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from os.path import splitext
import pandas as pd

base_path = r"/Users/Rydia/Dropbox (個人)/Lab/publications_conferences/SIEVE-Score_SciRep/"
auc_path = base_path + "data/auc/"
ef_path = base_path + "data/ef/"

target_list = open(base_path + "data/all.txt", "r").read().split("\n")
target_list = [x for x in target_list if len(x)>0]

aucs_result = []
aucs_glide = []
ef10s_result = []
ef10s_glide = []
ef1s_result = []
ef1s_glide = []


for target in target_list:
    auc_data = pd.read_csv(auc_path + target + ".csv", header=None)
    aucs_result.append(auc_data[0][5])
    aucs_glide.append(auc_data[0][6])

    ef_result = pd.read_csv(ef_path + target + "_result.txt", header=None)
    ef1s_result.append(ef_result[1][0])
    ef10s_result.append(ef_result[1][1])

    ef_glide = pd.read_csv(ef_path + target + "_glide.txt", header=None)
    ef1s_glide.append(ef_glide[1][0])
    ef10s_glide.append(ef_glide[1][1])

return_list = [("Target", "AUC_proposed", "AUC_Glide", "EF1_proposed", "EF1_Glide", "EF10_proposed", "EF10_Glide")]
return_list = return_list + \
              list(zip(target_list,
                       aucs_result, aucs_glide,
                       ef1s_result, ef1s_glide,
                       ef10s_result, ef10s_glide))

print(return_list)
return_df = pd.DataFrame(return_list)
return_df.to_csv(base_path + "data/result.csv", sep=",", index=False, header=False)