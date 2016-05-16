import numpy as np
import logging
logger = logging.getLogger(__name__)

def scoring_main(data,labels,n_clusters,outputfile,
                 plus,minus,propose,zeroneg,
                 score_type, score_params):

    score = []

    for k in range(n_clusters):
        class_members = labels == k
        x = scoring(k, data[class_members], outputfile, plus, minus,
                    propose, zeroneg, score_type, score_params)
        logger.info('Calc SIEVE-Score for cluster %d of %d.' 
                     % (k+1, n_clusters))
        logger.debug(x)
        score=score+list(x)

    score=sorted(score,key=lambda x:float(x[2]))[::-1][:propose] #.astype(str)

    np.savetxt(outputfile, score, fmt="%s", delimiter=",")
    logger.info('Saved SIEVE-Score.')


def scoring(cluster_number, data, outputfile, plus, minus, propose,
            zeroneg, score_type, score_params):

    title = data[:, 0]
    ishit = data[:, 1].astype('int')
    xyz = data[:, 3:].astype('float')

    hits = np.array([True if x > 0 else False for x in ishit])
    if zeroneg == False:
        nonhits = np.array([True if x < 0 else False for x in ishit])
    elif zeroneg == True:
        nonhits = np.array([True if x <= 0 else False for x in ishit])
    # <= for DUD-E testset, < for others

    hits_xyz = data[hits, 3:].astype('float')
    numhits = ishit[hits].shape[0]
    numnonhits = ishit[nonhits].shape[0]
    nonhits_xyz = data[nonhits, 3:].astype('float')
    known_xyz = np.r_[hits_xyz, nonhits_xyz]

    # calculate score
    score = np.zeros(data.shape[0])

    if numhits+numnonhits > 0:
        for i in range(data.shape[0]):
            for j in range(len(known_xyz)):
                if j < numhits:
                    multiplier = plus
                else:
                    multiplier = minus

                score[i] += score_func(xyz[i], known_xyz[j], 
                                       score_type, multiplier, *score_params)
    
        if ishit.shape[0] == 1:
            # single-hit cluster
            logger.warning('cluster %d: is single-hit cluster.'
                           %cluster_number)

    else:
        #contains only unknown data
        logger.warning('cluster %d: No known data.'%cluster_number)


    #exclude self-comparison
    for i in range(len(score)):
        if ishit[i] > 0:
            score[i] -= score_func(xyz[i], xyz[i],
                                   score_type, plus, *score_params)
        elif ishit[i] < 0 or (zeroneg == True and ishit[i] == 0):
            score[i] -= score_func(xyz[i], xyz[i],
                                   score_type, minus, *score_params)


    saved = (np.argsort(score)[::-1])[:propose]  # reversed,cut
    result = np.dstack((saved, title[saved],
                        score[saved], ishit[saved]))[0]
    return result

def score_func(a, b, type, *params):

    dist = lambda x, y: np.linalg.norm(x-y, 2)
    try:
        if type == 'normal':
            multiplier, cutoff, dim = params
            return multiplier / (max(dist(a,b), cutoff)**dim)

        elif type == 'exp':
            multiplier, cutoff, scale = params
            return multiplier / (max(np.exp(dist(a,b)/scale), cutoff))

        elif type == 'Gaussian':
            multiplier, sigma = params
            return (multiplier / ((2 * np.pi * sigma**2)**0.5)
                * np.exp(-dist(a,b)**2/(2 * sigma**2)))
        else:
            logger.error("Invalid score_type in score_func!")
            quit()

    except ValueError:
        logger.error("Invalid Argument in Score func!")
        quit()

def set_score_params(x):
    type = x['score_type']
    if type == 'normal':
        return x['cutoff'], x['score_dim']

    elif type == 'exp':
        return x['cutoff'], x['scale']

    elif type == 'Gaussian':
        return x['sigma']
 
    else:
        logger.error("Invalid score_type in set_score_params!")
        quit()
