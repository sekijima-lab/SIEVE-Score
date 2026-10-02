import logging
import numpy as np
import pandas as pd

def multi_importance(args):
    logger = logging.getLogger(__name__)

    from interaction_io import load_interaction
    from models import random_forest
    cpdname, label, interaction_name, interactions = load_interaction(
        args.input, args.hits, args.ignore)
    if args.mode == "interaction":
        return

    logger.info('Read interaction data.')


    if args.use_docking_score:
        X = interactions
    else:
        X = interactions[:, :-1]
    y = np.array([1 if x > 0 else 0 for x in label])


    importances = []
    for i in range(args.n_iter):
        seed = None if args.random_state is None else args.random_state + i
        clf = random_forest(seed, args.nprocs)
        importances.append(clf.fit(X, y).feature_importances_)
        print(i)
    importances = pd.DataFrame(importances, columns=interaction_name if args.use_docking_score else interaction_name[:-1], index=None)
    importances.to_csv("multiple_importance.csv", sep=",", index=None)

    print('\n*****Process Complete.*****\n')
    logger.info('\n*****Process Complete.*****\n')


if __name__ == '__main__':
    # options.py
    from options import Input_func as Input

    args = Input()

    logger = logging.getLogger(__name__)
    logger.info("options:\n" + str(args))

    multi_importance(args)
