from sklearn.metrics import roc_curve, roc_auc_score
import csv
import matplotlib
import math
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_ROC(n_actives, n_decoys, title, fname, results,
             legend, show):

    plt.figure(figsize=(6, 6), dpi=150)
    plt.xlabel("False Positive rate", fontsize=14)
    plt.ylabel("True Positive rate", fontsize=14)
    plt.title("ROC Curve (" + title + ")", fontsize=14)

    results_file = results[0::3]
    results_type = results[1::3]
    results_name = results[2::3]

    ys = []
    scores = []
    data_types = []

    for i in range(len(results_file)):
        y, score, t = GetYScoreFromResult(results_file[i], results_type[i])
        ys.append(y)
        scores.append(score)
        data_types.append(t)

    aucs = []
    ef10s = []
    ef1s = []
    for i in range(len(results_file)):
        auc, ef10, ef1 = AddROCCurve(plt, ys[i], scores[i], i,
                                     results_name[i], n_actives,
                                     n_decoys, data_types[i])
        aucs.append(auc)
        ef10s.append(ef10)
        ef1s.append(ef1)

    SaveROCCurvePlot(plt, fname, show, legend, True)
    SaveAUCEF(fname, results, aucs, ef10s, ef1s)


def SaveAUCEF(fname, results, aucs, ef10s, ef1s):
    from os.path import splitext

    results_file = results[0::3]
    results_type = results[1::3]
    results_name = results[2::3]

    glide_auc = proposed_auc = 0
    with open(splitext(fname)[0] + "_result.txt", "w") as f_result:
        f_result.write(str(title) + "\n")
        for i in range(len(results_file)):
            f_result.write(results_name[i] + "_AUC, " + str(aucs[i]) + "\n")
            f_result.write(results_name[i] + "_EF10, " + str(ef10s[i]) + "\n")
            f_result.write(results_name[i] + "_EF1, " + str(ef1s[i]) + "\n")

            if "Glide" in results_name[i]:
                glide_auc = aucs[i]

            elif "proposed" in results_name[i]:
                proposed_auc = aucs[i]

        if proposed_auc == 0 or glide_auc == 0:
            pass
        elif proposed_auc - glide_auc > 0.001:
            f_result.write("proposed, win\n")
        elif proposed_auc - glide_auc < -0.001:
            f_result.write("proposed, lose\n")
        else:
            f_result.write("proposed, draw\n")


def GetYScoreFromResult(filename, datatype):
    y = []
    score = []

    data = csv.reader(open(filename, 'rb'), delimiter=',', quotechar='#')

    if ('dock' in datatype or
            'fp' in datatype):
        # print(filename)
        for line in data:
            try:
                score.append(float(line[1]))
            except ValueError:
                continue
                # legend line?

            if "ZINC" in line[0] or "TEMP" in line[0]:
                y.append(0)
            else:
                y.append(1)

    elif datatype == 'stat':
        for line in data:
            try:
                y.append(int(line[3]))
                score.append(float(line[2]))
            except ValueError:
                pass

    else:
        print('unknown datatype.')
        sys.exit()

    reverse = bool(datatype != 'dock')
    res = zip(y, score)
    res = sorted(res, key=lambda x: x[1], reverse=reverse)
    y = [x[0] for x in res]
    score = [x[1] for x in res]

    return y, score, datatype


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


def AddROCCurve(plt, actives, scores, i, label, n_actives, n_decoys, type):

    tpr, fpr = GetRates(actives, scores, n_actives, n_decoys)

    colors = "rgbcmyk"
    color = colors[i % len(colors)]
    # colors= ['b', 'm', 'r']

    if "dock" in type:
        scores = [-x for x in scores]  # reverse order

    auc_tmp = roc_auc_score(actives, scores)
    auc = round((auc_tmp * tpr[-2] * fpr[-2] +
                 (tpr[-2] + 1) * (1 - fpr[-2]) / 2.0), 3)

    try:
        ef_10 = tpr[:-1][int(math.ceil((n_actives + n_decoys) / 10))] * 10
    except IndexError:
        ef_10 = tpr[-2] * 10

    try:
        ef_1 = tpr[:-1][int(math.ceil((n_actives + n_decoys) / 100))] * 100
    except IndexError:
        ef_1 = tpr[-2] * 100

    label = label + ", auc=" + str(auc)

    if "dock" in type:
        linestyle = '--'
    elif i < len(colors):
        linestyle = '-'
    else:
        linestyle = ':'

    plt.plot(fpr, tpr, linestyle=linestyle,
             color=color, linewidth=2, label=label)

    return auc, ef_10, ef_1


def SaveROCCurvePlot(plt, fname, show, legend, randomline=True):

    if randomline:
        x = [0.0, 1.0]
        plt.plot(x, x, linestyle='--', color='k',
                 linewidth=2, label='random')

    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    if legend:
        plt.legend(fontsize=12, loc='best')

    plt.tight_layout()
    plt.savefig(fname)
    if show:
        plt.show()


def num_molecule(x):
    if x.isdigit():
        return int(x)
    else:
        try:
            from schrodinger import structure
        except ImportError:
            logger.exception("error in counting molecule. " +
                             "if you want to count molecules of file, " +
                             "please execute by $SCHRODINGER/run python.",
                             exc_info=True)
            quit()

        st = structure.StructureReader(x)
        return sum(1 for _ in st)


if __name__ == '__main__':
    import sys
    x = sys.argv
    if len(x) <= 5 and len(x) % 3 != 2:
        print("usage: n_actives, n_decoys, graph_title, graph_filename,"
              "legend, show, results(file, type, name)...")
    n_actives = num_molecule(x[1])
    n_decoys = num_molecule(x[2])
    title = x[3]
    fname = x[4]
    legend = x[5]
    if legend in ["False", "false", "0", "No", "no"]:
        legend = False
    else:
        legend = True

    show = x[6]
    if show in ["False", "false", "0", "No", "no"]:
        show = False
    else:
        show = True

    results = x[7:]

    plot_ROC(n_actives, n_decoys, title, fname,
             results, legend, show)

