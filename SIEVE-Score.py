import logging
import numpy as np

def SIEVE(x):
    from sklearn.cluster import AffinityPropagation
    from read_interaction import read_interaction
    from scoring import scoring_main, set_score_params
    from sklearn.decomposition import PCA
    import matplotlib.pyplot as plt

    inter_array = np.array(read_interaction(x['i'],x['hits']))

    data = inter_array[1:, :]
    interactions = np.array(data[:, 2:], dtype='float')
    af = AffinityPropagation().fit(interactions)
    cluster_centers_indices = af.cluster_centers_indices_
    n_clusters = len(af.cluster_centers_indices_)
    labels = af.labels_

    logging.info('Read interaction data.')
    logging.debug('interaction_array:\n'+str(inter_array))
    logging.info('num_clusters: '+str(n_clusters))

    #setting score parameters
    score_type = x['score_type']
    score_params = set_score_params(x)

    logging.info('Score parameters: '+score_type+', '+str(score_params))

    _ = scoring_main(data,labels,n_clusters,x['o'],x['p'],x['m'],
                     x['propose'],x['zeroneg'],
                     score_type, score_params)

    pca = PCA(n_components=2)
    X = pca.fit(interactions).transform(interactions)
    from itertools import cycle
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for k, col in zip(range(n_clusters), colors):
        class_members = labels == k
        cluster_center = X[cluster_centers_indices[k]]
        plt.plot(X[class_members, 0], X[class_members, 1], col + '.')
    
    plt.title('Estimated number of clusters: %d' % n_clusters)
    plt.show()

    print('\n*****Process Complete.*****\n')
    logging.info('\n*****Process Complete.*****\n')

if __name__ == '__main__':
    import sys
    from options import Input_func as Input

    # options.py
    default_option = {'i': [], 'o': None, 'hits': None,
                      'cutoff': 1,'zeroneg': False,
                      'p': 1, 'm': -1,'p_opt': False,'opt_coef': 1,
                      'actives': None, 'decoys': None, 'propose': 100000,
                      'score_dim': 1.0, 'log': None, 'score_type':'normal'}

    if "-file" in sys.argv:
        option_file = sys.argv[sys.argv.index("-file") + 1]
        o = open(option_file, 'r').readlines()
        o = [s.replace("\n", "").strip() for s in o]
        o = o + sys.argv
    else:
        o = sys.argv
    
    options = Input(default_option,o)
    
    logger = logging.getLogger(__name__)
    logging.info("\noptions:"+str(options))

    # main function
    SIEVE(options)
