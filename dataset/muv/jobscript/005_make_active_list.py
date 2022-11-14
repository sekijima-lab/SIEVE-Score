#!/usr/bin/env python3

def muv_make_actives(filename, outputname=None):
    actives = []
    data = open(filename, "r").readlines()
    for l in data:
        smiles, cpdname = l.split(" ", 1)
        if cpdname.startswith("1.0 "):
            active = cpdname.strip("\n")+",1\n"
            actives.append(active)
    
    if outputname is None:
        from os import splitext
        outputname = splitext(filename)[0]+"_actives.txt"
    output = open(outputname, "w")
    output.writelines(actives)
    output.close()


if __name__=="__main__":
    proteins = ["CAG", "FAK", "FXI", "HIV",
                "HSP", "PKA", "SF1"]
    import os 
    
    for prot in proteins:
        os.chdir("../"+prot)
        muv_make_actives(prot+".smi", "actives.txt")

