#!/usr/bin/env python3

import logging
import numpy as np
import pandas as pd
import scoring

def load_interaction(f_name, hits, ignore=None):
    import os.path
    input_interaction = os.path.splitext(f_name)[0] + ".interaction"
    if os.path.exists(input_interaction):
        inter_array = pd.read_csv(input_interaction, index_col=0)
    else:
        from read_interaction import read_interaction
        data = np.array(read_interaction(f_name, hits))
        inter_array = pd.DataFrame(data[1:, 1:], index=data[1:,0].astype(str), columns=data[0,1:].astype(str), dtype=np.float64)

    # Split data
    
    interaction_name = inter_array.iloc[0, 2:].index.tolist()
    data = inter_array.iloc[1:, :]

    if ignore is not None:
        ignored = open(ignore, "r").read().split()
        data = data[~data.index.isin(ignored)]

    cpdname = data.index.values
    label = data.loc[:, "ishit"].values
    interactions = data.iloc[:, 1:].values
        
    return cpdname, label, interaction_name, interactions

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
