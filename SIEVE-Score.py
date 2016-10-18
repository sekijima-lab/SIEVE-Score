import logging
import numpy as np


def SIEVE(args):
    from sklearn.cluster import KMeans
    from read_interaction import read_interaction

    logger = logging.getLogger(__name__)

    inter_array = np.array(read_interaction(args.input, args.hits))

    data = inter_array[1:, :]
    interactions = np.array(data[:, 2:], dtype='float')

    logger.info('Read interaction data.')
    logger.debug('interaction_array:\n' + str(inter_array))

    # setting score parameters, eval
    from scoring import scoring_eval
    cpdname, score, label = scoring_eval(data, args)

    from plotting import plotting
    plotting(args, data, cpdname, score, label)
    print('\n*****Process Complete.*****\n')
    logger.info('\n*****Process Complete.*****\n')

if __name__ == '__main__':

    # options.py
    from options import Input_func as Input
    args = Input()

    logger = logging.getLogger(__name__)
    logger.info("options:\n" + str(args))

    SIEVE(args)
