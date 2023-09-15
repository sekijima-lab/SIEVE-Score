import csv
import sys
from os.path import splitext

import numpy as np


def test():
    """For auto-testing"""
    n_actives = 100
    n_decoys = 4900
    
    y = [1]*100+[0]*4900
    print(calc_ef(y, n_actives, n_decoys, threshold=0.1)) # 10
    print(calc_ef(y, n_actives, n_decoys, threshold=0.01)) # 50

    y = [0]*4900+[1]*100
    print(calc_ef(y, n_actives, n_decoys, threshold=0.1)) # 0
    print(calc_ef(y, n_actives, n_decoys, threshold=0.01)) # 0

    y = [1]*50+[0]*4900+[1]*50
    print(calc_ef(y, n_actives, n_decoys, threshold=0.1)) # 5
    print(calc_ef(y, n_actives, n_decoys, threshold=0.01)) # 50 


def main(n_actives, n_decoys, f_active, f_result, f_out=None):

    if "mae" in splitext(f_result)[1]:
        y, score = get_y_score_from_glide(f_result, f_active)
    else:
        y, score = get_y_score_from_result(f_result, f_active)

    # sort y, score
    print(y, score)

    ef10 = calc_ef(y, n_actives, n_decoys, 0.1)
    ef1 = calc_ef(y, n_actives, n_decoys, 0.01)

    if f_out != sys.stdout:
        print("EF_01," + str(ef1))
        print("EF_10," + str(ef10))

    with open(f_out, 'w') as f_out:
        f_out.write("EF_01," + str(ef1) + "\n")
        f_out.write("EF_10," + str(ef10) + "\n")


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

    # reverse
    y = np.array(ys)
    score = np.array(scores)
    y = y[np.argsort(score)][::-1]
    score = np.sort(score)[::-1]

    return y, score


def get_y_score_from_glide(f_result, f_active):

    try:
        from schrodinger import structure
    except ImportError:
        print("if you want to use this function, " +
              "please execute from $SCHRODINGER/run python.")
        quit()

    reader = structure.StructureReader(f_result)
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

    # not reverse
    y = np.array(ys)
    score = np.array(scores)
    y = y[np.argsort(score)]
    score = np.sort(score)

    return y, score


def count_actives_decoys(f_glide, f_active):
    try:
        from schrodinger import structure
    except ImportError:
        print("if you want to count molecules of glide result, " +
              "please execute from $SCHRODINGER/run python.")
        quit()

    reader = structure.StructureReader(f_glide)
    actives = get_active_from_activefile(f_active)
    n_actives = n_decoys = 0

    for st in reader:
        prop = st.property

        if 'r_i_docking_score' not in prop.keys():  # protein, or error?
            continue
        elif prop['s_m_title'] in actives:
            n_actives += 1
        else:
            n_decoys += 1

    reader.close()
    return n_actives, n_decoys


def calc_ef(sorted_y, n_actives=None, n_decoys=None, threshold=0.1):
    """Calculate enrichment factor.

    Args:
        sorted_y (array or list-like): Your screen result Y. Higher score should be first.
        n_actives ([type], optional): Number of actives. If not supplied, automatically use number of label==1.
        n_decoys ([type], optional): Number of actives. If not supplied, automatically use number of label==0.
        threshold (float, optional): x% enrichment factor in 0--1 value. Defaults to 0.1.

    Returns:
        float: (threshold*100)% enrichment factor.
    """
    if not 0 <= threshold <= 1:
        print("Error in calc_ef: threshold must be 0<=x<=1")
        quit()

    y = np.array(sorted_y)
    if n_actives is None:
        n_actives = np.sum(y==1)
    if n_decoys is None:
        n_decoys = np.sum(y==0)
    total = n_actives + n_decoys

    random_rate = float(n_actives) / total
    screen_range = int(np.ceil(total * threshold))
    screen_rate = float(sum(y[:screen_range])) / screen_range

    return screen_rate / random_rate


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(prog="SIEVE-Score_v1.5-analysis",
                                     description="Calculate EF from result",
                                     fromfile_prefix_chars='@')
    parser.add_argument("input", help="input Glide pv file or result csv")
    parser.add_argument("-a", "--active", required=True,
                        help="active definition file")
    parser.add_argument("-o", "--output", nargs="?",
                        default=sys.stdout, help="output file")
    parser.add_argument("-n", "--number", nargs=2, default=None,
                        help="number of actives, decoys (2 integers), prior to -d")
    parser.add_argument("-d", "--dock", nargs="?",
                        default=None, help="dock result file for number of actives, decoys")

    args = parser.parse_args(sys.argv[1:])

    if args.number is not None:
        n_actives, n_decoys = args.number
        n_actives = int(n_actives)
        n_decoys = int(n_decoys)
    elif args.dock is not None:
        n_actives, n_decoys = count_actives_decoys(args.dock, args.active)
    else:
        print("Either of -n or -d is needed.")
        quit()

    return args, n_actives, n_decoys


def init():
    args, n_actives, n_decoys = parse_args()

    main(n_actives, n_decoys, args.active, args.input, args.output)


if __name__ == '__main__':
    init()
