from schrodinger import structure
import numpy as np
from os.path import splitext
import logging

logger = logging.getLogger(__name__)


def read_interaction_file(inputfile, residues, hits, result):
    reader = structure.StructureReader(inputfile)
    for st in reader:
        prop = st.property
        inter = []

        if 'r_i_docking_score' not in prop.keys():  # protein? error?
            continue

        inter.append(prop['s_m_title'])
        if hits is not None and prop['s_m_title'] in hits.index:
            i = hits.loc[prop["s_m_title"], :].iloc[0]
            inter.append(i)
        else:
            inter.append(0)

        for x in residues:
            if x == "docking_score":
                inter.append(float(prop["r_i_docking_score"]))
            else:
                vdW = float(prop['r_glide_res:' + x + '_vdw'])
                coul = float(prop['r_glide_res:' + x + '_coul'])
                hbond = float(prop['r_glide_res:' + x + '_hbond'])
                inter = inter + [vdW, coul, hbond]

        result.append(inter)

    reader.close()
    return result


def read_label_mae(hitsfile, label):
    reader = structure.StructureReader(hitsfile)
    hits = []
    for st in reader:
        prop = st.property
        if 'r_i_docking_score' in prop.keys():  # protein? error?
            hits.append((prop['s_m_title'], label))

    dt = np.dtype([('title', np.str_, 64), ('ishit', np.int64, 1)])
    hits = np.array(hits, dtype=dt)
    if hits.shape == ():
        hits = np.array([hits])

    return hits


def read_label_smi(hitsfile, label):
    reader = open(hitsfile).readlines()
    hits = []
    for st in reader:
        molname = st.replace("\n","").split(" ")[-1]
        hits.append((molname, label))
    dt = np.dtype([('title', np.str_, 64), ('ishit', np.int64, 1)])
    hits = np.array(hits, dtype=dt)
    if hits.shape == ():
        hits = np.array([hits])
    return hits

def read_interaction(inputfile, hitsfile):
    reader = structure.StructureReader(inputfile)
    st = next(reader)
    while 'r_i_docking_score' not in st.property.keys():
        try:
            st = next(reader)  # skip protein
        except StopIteration:
            logger.exception("error in reading maegz file. " +
                             "maybe it does not contain interaction data.",
                             exc_info=True)
            quit()

    residues = []
    result = []
    legend = ['# title', 'ishit']

    for x in st.property.keys():
        if x.startswith('r_glide_res:') and\
           x.endswith('_Eint'):
            res = x.replace('r_glide_res:', '').replace('_Eint', '')
            residues.append(res)
    residues.append("docking_score")
    residues = sorted(residues)

    for res in residues:
        if res != "docking_score":
            res_str = str(res)
            legend = legend + [res_str+"_vdw", res_str+"_coul", res_str+"_hbond"]
        else:
            legend.append("docking_score")

    result.append(legend)

    reader.close()
    
    if hitsfile is None:
        hits = None
    elif splitext(hitsfile)[1] in [".mae", ".maegz", ".gz"]:
        hits = read_label_mae(hitsfile, 1)
    elif splitext(hitsfile)[1] in [".ism", ".smi"]:
        hits = read_label_smi(hitsfile, 1)
        # print(hits)
    else:
        # print(hitsfile)
        try:
            import pandas as pd
            hits = pd.read_csv(hitsfile, sep=",", comment="#",
                               header=None, index_col=0,
                               names=["ishit"])
            if hits.shape == ():
                hits = np.array([hits])
        except IOError:
            hits = None

    result = read_interaction_file(inputfile, residues, hits, result)
    # print(hits)

    int_file = splitext(inputfile)[0]+".interaction"
    np.savetxt(int_file, result, delimiter=",", fmt="%s")

    return result
