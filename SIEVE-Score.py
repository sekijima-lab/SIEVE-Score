import logging
import numpy as np

def SIEVE(args):
    from sklearn.cluster import KMeans
    from read_interaction import read_interaction

    logger = logging.getLogger(__name__)

    inter_array = np.array(read_interaction(args.input,args.hits))

    data = inter_array[1:, :]
    interactions = np.array(data[:, 2:], dtype='float')

    n_clusters = args.num_cluster
    km = KMeans(n_clusters=n_clusters,
                max_iter=1000,n_jobs=-1).fit(interactions)
    labels = km.labels_

    logger.info('Read interaction data.')
    logger.debug('interaction_array:\n'+str(inter_array))
    logger.info('num_clusters: '+str(n_clusters))

    #setting score parameters
    from scoring import scoring_main, set_score_params
    score_type = args.score_type
    score_params = set_score_params(args)

    logger.info('Score parameters: '+score_type+', '+str(score_params))

    score = scoring_main(data,labels,n_clusters,args,
                     score_type, score_params)

    from plotting import plotting
    plotting(args, data, km, score)
    print('\n*****Process Complete.*****\n')
    logger.info('\n*****Process Complete.*****\n')

if __name__ == '__main__':

    # options.py
    from options import Input_func as Input
    args = Input()
    
    logger = logging.getLogger(__name__)
    logger.info("options:\n"+str(args))


    SIEVE(args)
