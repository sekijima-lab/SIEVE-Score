import numpy as np
import logging
from multiprocessing import Pool

logger = logging.getLogger(__name__)

def argwrapper(x):
    return x[0](*x[1:])

def dist_euclidian(a, b):
    return np.linalg.norm(a-b, 2)

def dist_mahalanobis(a, b):
    import sklearn.neighbors.DistanceMetric
    dist = DistanceMetric.get_metric('mahalanobis')
    return 

def scoring_cluster(n_clusters, cluster_number, data, labels, 
                    args, score_func):    
    class_members = labels == cluster_number
    data_in_cluster = data[class_members]

    result = scoring(cluster_number, data_in_cluster, args, score_func)
    logger.info('Calc SIEVE-Score for cluster %d of %d.'
                % (cluster_number + 1, n_clusters))
    logger.debug(result)

    return result


def scoring(cluster_number, data, args, score_func):

    plus = args.plus
    minus = args.minus
    propose = args.propose
    zeroneg = args.zeroneg

    title = data[:, 0]
    ishit = data[:, 1].astype('int')
    xyz = data[:, 3:].astype('float')

    hits = np.array([True if x > 0 else False for x in ishit])
    if zeroneg is False:
        nonhits = np.array([True if x < 0 else False for x in ishit])
    elif zeroneg is True:
        nonhits = np.array([True if x <= 0 else False for x in ishit])
    # <= for DUD-E testset, < for others

    hits_xyz = data[hits, 3:].astype('float')
    numhits = ishit[hits].shape[0]
    numnonhits = ishit[nonhits].shape[0]
    nonhits_xyz = data[nonhits, 3:].astype('float')
    known_xyz = np.r_[hits_xyz, nonhits_xyz]

    # calculate score
    score = np.zeros(data.shape[0])

    if numhits + numnonhits > 0:
        for i in range(data.shape[0]):
            for j in range(len(known_xyz)):
                if j < numhits:
                    multiplier = plus
                else:
                    multiplier = minus

                score[i] += score_func.calculate(multiplier,
                                                 xyz[i], known_xyz[j])

        if ishit.shape[0] == 1:
            # single-hit cluster
            logger.warning('cluster %d: is single-hit cluster.'
                           % cluster_number)

    else:
        # contains only unknown data
        logger.warning('cluster %d: No known data.' % cluster_number)

    # exclude self-comparison
    for i in range(len(score)):
        if ishit[i] > 0:
            score[i] -= score_func.calculate(multiplier,
                                             xyz[i], known_xyz[j])
        elif ishit[i] < 0 or (zeroneg is True and ishit[i] == 0):
            score[i] -= score_func.calculate(multiplier,
                                             xyz[i], known_xyz[j])

    saved = (np.argsort(score)[::-1])[:propose]  # reversed,cut
    result = np.dstack((saved, title[saved],
                        score[saved], ishit[saved]))[0]
    return result


def save_to_maegz(inputfile, labels, score):
    import schrodinger.structure as structure

    reader = structure.StructureReader(inputfile)
    write_maegz = inputfile.rstrip(".maegz") + "_result.maegz"
    writer = structure.StructureWriter(write_maegz)

    index = 0
    for st in reader:
        prop = st.property

        if 'r_i_docking_score' not in prop.keys():  # protein? error?
            writer.append(st)
            continue
        else:
            st.property['i_Clustering_Cluster#'] = labels[index]
            st.property['i_Clustering_PCAScore'] = score[index]
            writer.append(st)
        index += 1

    reader.close()
    writer.close()
    print('Saved Cluster info to ' + write_maegz)


def scoring_main(data, labels, n_clusters, args):

    clustering = args.clustering
    outputfile = args.output
    score_type = args.score_type

    score = []

    from scoring_function import set_score_func

    #CHANGE HERE
    score_func = set_score_func(args, dist)
    logger.info(str(score_func))

    if clustering is True:
        pool = Pool(processes=4)

        score_nest = pool.map(argwrapper,
                              [(scoring_cluster, n_clusters, i,
                                data, labels, args, score_func)
                               for i in range(n_clusters)])

        score = [_ for inner in score_nest for _ in inner]
        # flatten list

    else:
        # clustering is False
        x = scoring(0, data, args, score_func)
        logger.info('Calc SIEVE-Score for cluster %d of %d.'
                    % (1, n_clusters))
        logger.debug(x)
        score = score + list(x)

    score = sorted(score, key=lambda x: float(x[2]))[::-1][:args.propose]

    np.savetxt(outputfile, score, fmt="%s", delimiter=",")
    logger.info('Saved SIEVE-Score.')

    return score
