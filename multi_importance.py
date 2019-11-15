import logging
import numpy as np
import pandas as pd

def multi_importance(args):
    logger = logging.getLogger(__name__)

    # read interaction
    import os.path
    input_interaction = os.path.splitext(args.input[0])[0] + ".interaction"
    if len(args.input) == 1 and os.path.exists(input_interaction):
        inter_array = np.genfromtxt(input_interaction, comments=None, delimiter=",", dtype=None)
        inter_array = inter_array.tolist()
        inter_array = np.asarray(inter_array)
    else:
        from read_interaction import read_interaction
        inter_array = np.array(read_interaction(args.input, args.hits))

    if args.mode == "interaction":
        logger.info("Saved interaction data.")
        logger.info("\n****Process Complete.****")
        quit()

    # Split data
    interaction_name = inter_array[0, 2:].astype('str')
    data = inter_array[1:, :]

    cpdname = data[:, 0].astype('str')
    label = data[:, 1].astype('float')
    interactions = data[:, 2:].astype('float')

    logger.info('Read interaction data.')
    logger.debug('interaction_array:\n' + str(inter_array))

    if args.use_docking_score:
        X = interactions
    else:
        X = interactions[:, :-1]
    y = np.array([1 if x > 0 else 0 for x in label])
    from sklearn.ensemble import RandomForestClassifier as RFC
    clf = RFC(n_estimators=1000, criterion='gini', max_features=6, n_jobs=args.nprocs)
    model_names = 'SIEVE-Score ' + args.model
    
    importances = []
    for i in range(args.n_iter):
        importances.append(clf.fit(X, y).feature_importances_)
        print(i)
    importances = pd.DataFrame(importances, columns=interaction_name[:-1], index=None)
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
