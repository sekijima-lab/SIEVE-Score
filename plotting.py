import logging
import numpy as np
from sklearn.decomposition import PCA
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def plotting(args, data, cpd_name, score, label):
    interactions = np.array(data[:, 2:], dtype='float')
    cpd_name = data[:, 0]

    pca = PCA(n_components=2)
    X = pca.fit(interactions).transform(interactions)
    # from itertools import cycle
    # colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')

    print(score, label)
    class_members = (label > 0)
    color_map = 'hot'
    plt.scatter(X[class_members, 0], X[class_members, 1],
                c=score[class_members], marker='o', cmap=color_map, vmin=0, vmax=1)

    class_members = (label <= 0)
    color_map = 'hot'
    plt.scatter(X[class_members, 0], X[class_members, 1],
                c=score[class_members], marker='x', cmap=color_map, vmin=0, vmax=1)

    plt.colorbar()

    plt.xlabel("PC1")
    plt.ylabel("PC2")

    # cpd_name = data[:,0]
    if args.annotate:
        # debug mapping, annotating
        for index, xy in enumerate(X[:, 0:2]):
            score_cpd = score[index]

            plt.annotate('%.3f, %s' % (score_cpd, cpd_name[index]),
                         xy=xy, textcoords='data', size='small')

    plt.axes().set_aspect("equal", "datalim")

    plt.title(args.title + " in PCA space")

    from os.path import splitext
    figname = splitext(args.output)[0] + '.png'
    plt.savefig(figname, dpi=300)
    logger.info('saved figure to %s' % figname)
