import logging
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
logger = logging.getLogger(__name__)


def plotting(args, data, cpdname, score, label):

    interactions = np.array(data[:, 2:], dtype='float')
    cpdname = data[:, 0]

    pca = PCA(n_components=2)
    X = pca.fit(interactions).transform(interactions)

    #from itertools import cycle
    #colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')

    class_members = (label > 0)
    plt.plot(X[class_members, 0], X[class_members, 1], 'b' + 'o')
    class_members = (label <= 0)
    plt.plot(X[class_members, 0], X[class_members, 1], 'g' + 'x')

    plt.xlabel("PC1")
    plt.ylabel("PC2")

    # cpdname = data[:,0]
    if args.annotate:
        # debug mapping, annotating
        for index, xy in enumerate(X[:, 0:2]):
            score_vector = [v for v in score if v[1] == cpdname[index]][0]
            if score_vector == []:
                # not found
                score_cpd = ""
            else:
                score_cpd = float(score_vector[2])

            plt.annotate('%.3f, %s' % (score_cpd, cpdname[index]),
                         xy=xy, textcoords='data', size='small')

    plt.axes().set_aspect("equal", "datalim")

    plt.title(args.title)

    from os.path import splitext
    figname = splitext(args.output)[0] + '.png'
    plt.savefig(figname, dpi=300)
    logger.info('saved figure to %s' % figname)
