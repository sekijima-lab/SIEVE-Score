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

        if hits is not None and prop['s_m_title'] in hits['title']:
            place = np.where(hits['title'] == prop['s_m_title'])
            i = int(place[0][0])
            inter.append(hits[i][1])
        else:
            inter.append(0)

        for x in residues:
            if x == "docking_score":
                inter.append(float(prop["r_i_docking_score"]))
            else:
                vdW = float(prop['r_glide_res:' + x + '_vdW'])
                coul = float(prop['r_glide_res:' + x + '_coul'])
                hbond = float(prop['r_glide_res:' + x + '_hbond'])
                inter = inter + [vdW, coul, hbond]

        result.append(inter)

    reader.close()
    return result


def read_label(hitsfile, label):
    reader = structure.StructureReader(hitsfile)
    for st in reader:
        hits = []
        prop = st.property
        if 'r_i_docking_score' in prop.keys():  # protein? error?
            hits.append((prop['s_m_title'], label))

    dt = np.dtype([('title', np.str_, 64), ('ishit', np.int64, 1)])
    hits = np.array(hits, dtype=dt)
    if hits.shape == ():
        hits = np.array([hits])

    return hits


def read_interaction(inputfiles, hitsfile):
    reader = structure.StructureReader(inputfiles[0])

    st = reader.next()
    while 'r_i_docking_score' not in st.property.keys():
        try:
            st = reader.next()  # skip protein
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

    #residues = sorted(residues)
    for res in residues:
        res_str = str(res)
        legend = legend + [res_str+"_vdW", res_str+"coul", res_str+"hbond"]

    result.append(legend)

    reader.close()

    if splitext(hitsfile)[1] in [".mae", ".maegz", ".gz"]:
        hits = read_label(hitsfile, 1)
        # print(hits)
    else:
        # print(hitsfile)
        try:
            hits = np.loadtxt(hitsfile, delimiter=',', comments='#',
                              dtype={'names': ('title', 'ishit'),
                                     'formats': ('S64', '<i2')})
            if hits.shape == ():
                hits = np.array([hits])
        except IOError:
            hits = None

    for f in inputfiles:
        result = read_interaction_file(f, residues, hits, result)
    # print(hits)
    return result
