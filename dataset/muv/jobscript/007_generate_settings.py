#!/usr/bin/env python3

def generate_settings(target_name):
    header=[
    "-i",
    "glide-dock_best_pv.maegz",
    "-l",
    "SIEVE-Score.log",
    "--hits",
    "actives.txt",
    "-z",
    "-d",
    "--n_iter",
    "5"]

    compounds_file = "../"+target_name+"/"+target_name+".smi"
    data = open(compounds_file, "r").read()
    actives = data.count("1.0")
    decoys = data.count("0.0")

    settings = ["--active", str(actives),
                "--decoy", str(decoys),
                "-t", "SIEVE-Score: "+target_name]

    output = open("../"+target_name+"/SIEVE-Score_settings.txt", "w")
    output.write("\n".join(header+settings))
    output.close()


if __name__=="__main__":
    proteins = ["CAG", "FAK", "FXI", "HIV",
                "HSP", "PKA", "SF1"]

    for prot in proteins:
        generate_settings(prot)
