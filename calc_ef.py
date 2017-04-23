from sklearn.metrics import roc_curve, roc_auc_score
import csv
from os.path import splitext
import sys
import pandas as pd
import numpy as np


def calc_EF(n_actives, n_decoys, title, fname, onlyAUC, results):
    results_file = results
    names = [splitext(x)[0] for x in results_file]

    aucs = []
    ef10s = []
    ef1s = []

    for i in range(len(results_file)):
        y, score = GetYScoreFromResult(results_file[i])
        auc_tmp = roc_auc_score(y, score)

        tpr, fpr = GetRates(y, score, n_actives, n_decoys)
        auc = round((auc_tmp * tpr[-2] * fpr[-2] +
                     (tpr[-2] + 1) * (1 - fpr[-2]) / 2.0), 3)

        ef10 = tpr[len(tpr) // 10] * 10
        ef1 = tpr[len(tpr) // 100] * 100

        aucs.append(auc)
        ef10s.append(ef10)
        ef1s.append(ef1)

    SaveAUCEF(fname, results, names, onlyAUC, aucs, ef10s, ef1s)


def GetYScoreFromResult(filename):
    res = []

    data = csv.reader(open(filename, 'rb'), delimiter=',', quotechar='#')

    for line in data:
        res.append((int(line[3]), float(line[2])))
        # y=line[3], score=line[2]

    res = sorted(res, key=lambda x: x[1], reverse=True)
    y = [x[0] for x in res]
    score = [x[1] for x in res]

    return y, score


def GetRates(y, scores, n_actives, n_decoys):

    tpr = [0.0]  # true positive rate
    fpr = [0.0]  # false positive rate

    foundactives = 0.0
    founddecoys = 0.0
    for idx, score in enumerate(scores):
        if y[idx] == 1:
            foundactives += 1.0
        else:
            founddecoys += 1.0

        tpr.append(foundactives / float(n_actives))
        fpr.append(founddecoys / float(n_decoys))

    tpr.append(1.0)
    fpr.append(1.0)  # add [1.0, 1.0]

    return tpr, fpr


def SaveAUCEF(fname, results, names, onlyAUC, aucs, ef10s, ef1s):

    with open(fname, "w") as f_result:
        f_result.write(str(title) + "\n")
        for i in range(len(results)):
            f_result.write(names[i] + "_AUC, " +
                           str(aucs[i]) + "\n")

            if not onlyAUC:
                f_result.write(names[i] + "_EF10, " + str(ef10s[i]) + "\n")
                f_result.write(names[i] + "_EF1, " + str(ef1s[i]) + "\n")

    x = [aucs[7 * i:7 * (i + 1)] for i in range(7)]

    df = pd.DataFrame(x)
    df.columns = [0.2, 0.5, 1, 2, 5, 10, 20]
    df.index = [0.125, 0.25, 0.5, 1, 2, 4, 8]
    df.to_csv(splitext(fname)[0] + ".csv")

if __name__ == '__main__':
    import sys
    x = sys.argv
    if len(x) <= 5 and len(x) % 3 != 2:
        print("usage: n_actives, n_decoys, graph_title, graph_filename,"
              "legend, show, glidefile, results(file, type, name)...")

    n_actives = int(x[1])
    n_decoys = int(x[2])
    title = x[3]
    fname = x[4]
    if x[5] in ["False", "false", "No", "no", "0", False]:
        onlyAUC = False
    else:
        onlyAUC = True

    results = x[6:]

    calc_AUC(n_actives, n_decoys, title, fname, onlyAUC, results)
