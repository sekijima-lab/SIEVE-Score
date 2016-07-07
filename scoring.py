import numpy as np
import logging
from multiprocessing import Pool

logger = logging.getLogger(__name__)

def argwrapper(x):
    return x[0](*x[1:])

def scoring_cluster(n_clusters, cluster_number, data, labels, args):
    class_members = labels == cluster_number
    data_in_cluster = data[class_members]

    from scoring_function import set_score_func
    #score_func = set_score_func(args, dist=YOUR_DISTANCE_FUNCTION)
    score_func = set_score_func(args)

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

    saved = np.argsort(score)[::-1]  # reversed,cut
    result = np.dstack((saved, title[saved],
                        score[saved], ishit[saved]))[0]
    return result


def scoring_main(data, labels, n_clusters, args):

    clustering = args.clustering
    outputfile = args.output
    score_type = args.score_type

    score = []

    if clustering is True:
        pool = Pool(processes=args.nprocs)

        score_nest = pool.map(argwrapper,
                              [(scoring_cluster, n_clusters, i,
                                data, labels, args)
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


    #score: [[id, title, score, label], ...]
    score = sort_delete_duplicate(score)
    score = score[:args.propose]

    np.savetxt(outputfile, score, fmt="%s", delimiter=",")
    logger.info('Saved SIEVE-Score.')

    return score


def sort_delete_duplicate(score):
    # sort by title, score
    score = sorted(score, key=lambda x: float(x[2]), reverse=True)
    score = sorted(score, key=lambda x: x[1], reverse=True)

    ret = []
    for d in score:
        cpdname = d[1]
        if len(ret)>0:
            lastcpdname = ret[-1][1]
        else:
            lastcpdname = None

        if cpdname != lastcpdname:
            ret.append(d)

    # re-sort by score
    ret = sorted(ret, key=lambda x: float(x[2]), reverse=True)
    return ret
