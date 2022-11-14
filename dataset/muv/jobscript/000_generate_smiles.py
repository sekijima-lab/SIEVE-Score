#!/usr/bin/env python3

import pandas as pd
import os

def writesmi(dataobj, rowname, outputname):
    row = dataobj.dropna(subset=[rowname]).loc[:, ["smiles", rowname, "mol_id"]]
    row.to_csv("../"+outputname+"/"+outputname+".smi", index=None, sep=" ", header=None)

#dataset: https://github.com/deepchem/deepchem/raw/master/datasets/muv.csv.gz

data = pd.read_csv("muv.csv", sep=",")
pdbs = ["MUV-548", "MUV-600", "MUV-652", "MUV-712", 
        "MUV-810", "MUV-832", "MUV-846"]
outputnames = ["PKA", "SF1", "HIV", "HSP",
               "FAK", "CAG", "FXI"]

for p, o in zip(pdbs, outputnames):
    try:
        os.mkdir("../"+o)
    except FileExistsError:
        pass
    writesmi(data, p, o)
