import logging
import numpy as np


def sieve(args):
    logger = logging.getLogger(__name__)

    # read interaction
    import os.path
    input_interaction = os.path.splitext(args.input[0])[0] + ".interaction"
    if len(args.input) == 1 and os.path.exists(input_interaction):
        inter_array = np.genfromtxt(input_interaction,
                                    delimiter=",", dtype=None)
        inter_array = inter_array.tolist()
        inter_array = np.asarray(inter_array)
    else:
        from read_interaction import read_interaction
        inter_array = np.array(read_interaction(args.input, args.hits))
        if args.inter_only:
            logger.info("Saved interaction data.")
            logger.info("\n****Process Complete.****")
            quit()

    # Split data
    # header = inter_array[0, :]
    data = inter_array[1:, :]

    cpdname = data[:, 0]
    label = data[:, 1].astype('float')
    interactions = data[:, 2:].astype('float')

    logger.info('Read interaction data.')
    logger.debug('interaction_array:\n' + str(inter_array))

    # Calc SIEVE-Score
    from scoring import scoring_eval, scoring_param_search, scoring_compareSVMRF
    # scoring_param_search(cpdname, label, interactions, args)
    scoring_eval(cpdname, label, interactions, args)
    # scoring_compareSVMRF(cpdname, label, interactions, args)

    # plot on PCA space
    # from plotting import plotting
    # plotting(args, data, cpdname, score, label)

    print('\n*****Process Complete.*****\n')
    logger.info('\n*****Process Complete.*****\n')


if __name__ == '__main__':
    # options.py
    from options import Input_func as Input

    args = Input()

    logger = logging.getLogger(__name__)
    logger.info("options:\n" + str(args))

    sieve(args)
