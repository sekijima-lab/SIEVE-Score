import csv
from os.path import splitext
import numpy as np


def main(n_actives, n_decoys, f_active, f_result, f_out=None):

    if "mae" in splitext(f_result)[1]:
        y, score = get_y_score_from_glide(f_result, f_active)
    else:
        y, score = get_y_score_from_result(f_result, f_active)

    # sort y, score
    y = np.array(y)
    score = np.array(score)
    y = y[np.argsort(score)][::-1]
    score = np.sort(score)[::-1]
    print(y, score)

    ef10 = calc_ef(y, n_actives, n_decoys, 0.1)
    ef1 = calc_ef(y, n_actives, n_decoys, 0.01)

    print("EF_01, "+str(ef1))
    print("EF_10, "+str(ef10))

    if f_out is not None:
        with open(f_out, "w") as f:
            f.write("EF_01," + str(ef1) + "\n")
            f.write("EF_10," + str(ef10) + "\n")


def get_active_from_activefile(f_active):

    actives = csv.reader(open(f_active, 'rb'), delimiter=',', quotechar='#')
    ret = []

    for line in actives:
        if float(line[1]) > 0:
            ret.append(line[0])

    return ret


def get_y_score_from_result(f_result, f_active):

    data = csv.reader(open(f_result, 'rb'), delimiter=',', quotechar='#')
    actives = get_active_from_activefile(f_active)

    ys = []
    scores = []

    for line in data:
        name = line[0]
        val = float(line[1])

        if name in actives:
            ys.append(1)
        else:
            ys.append(0)
        scores.append(val)

    return ys, scores


def get_y_score_from_glide(f_result, f_active):

    from schrodinger import structure as sts

    reader = sts.StructureReader(f_result)
    actives = get_active_from_activefile(f_active)

    ys = []
    scores = []
    for st in reader:
        prop = st.property

        if 'r_i_docking_score' not in prop.keys():  # protein, or error?
            continue

        name = prop['s_m_title']
        val = float(prop["r_i_docking_score"])

        if name in actives:
            ys.append(1)
        else:
            ys.append(0)
        scores.append(val)

    reader.close()
    return ys, scores


def calc_ef(y, n_actives, n_decoys, threshold=0.1):

    if not 0 <= threshold <= 1:
        print("Error in calc_ef: threshold must be 0<x<1")
        quit()

    total = n_actives+n_decoys
    random_rate = float(n_actives) / total
    screen_range = int(np.ceil(total * threshold))

    screen_rate = float(sum(y[:screen_range])) / screen_range

    return screen_rate / random_rate


def init():
    import sys
    x = sys.argv
    if len(x) <= 4:
        print("usage: n_actives, n_decoys, f_active, f_result [, f_out]")
        quit()

    n_actives = int(x[1])
    n_decoys = int(x[2])
    f_active = x[3]
    f_result = x[4]

    if len(x) == 6:
        f_out = x[5]
    else:
        f_out = None

    main(n_actives, n_decoys, f_active, f_result, f_out)


if __name__ == '__main__':
    init()
