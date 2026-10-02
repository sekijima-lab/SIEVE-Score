#!/usr/bin/env python3

import logging
import scoring
from interaction_io import load_interaction

def sieve(args):
    logger = logging.getLogger(__name__)

    # read interaction
    cpdname, label, interaction_name, interactions = load_interaction(args.input, args.hits, args.ignore)

    if args.mode == "interaction":
        logger.info("Saved interaction data.")
        logger.info("\n****Process Complete.****")
        quit()
    elif args.mode == "screen":
        cpdname_t, label_t, interaction_name_t, interactions_t = load_interaction(args.testdata, args.testhits)

    logger.info('Read interaction data.')

    # Calc SIEVE-Score
    if args.mode == "paramsearch":
        logger.info("SIEVE: main: Do parameter search")
        scoring.scoring_param_search(cpdname, label, interaction_name, interactions, args)
    elif args.mode == "comparesvm":
        logger.info("SIEVE: main: Do compare SVM and RF")
        scoring.scoring_compareSVMRF(cpdname, label, interaction_name, interactions, args)
    elif args.mode == "cv":
        logger.info("SIEVE: main: Do evaluation")
        scoring.scoring_eval(cpdname, label, interaction_name, interactions, args)
    elif args.mode == "datasize":
        logger.info("SIEVE: main: Do datasize analysis")
        scoring.lessdata_auc_ef(cpdname, label, interaction_name, interactions, args)
    elif args.mode == "screen":
        from screening import screening
        logger.info("SIEVE: main: Do screening")
        screening(cpdname, label, interaction_name, interactions,
                  cpdname_t, label_t, interaction_name_t, interactions_t, args)
    else:
        logger.info("SIEVE: main: Unknown mode?")
        quit()

    print('\n*****Process Complete.*****\n')
    logger.info('\n*****Process Complete.*****\n')


if __name__ == '__main__':
    # options.py
    from options import Input_func as Input

    args = Input()

    logger = logging.getLogger(__name__)
    logger.info("options:\n" + str(args))

    sieve(args)
