import numpy as np
import pandas as pd
import logging
from os.path import splitext
from models import random_forest, support_vector

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)
cv = 5  # define k of k-fold CV


def lessdata_auc_ef(cpd_names, label_data, interaction_name, features, args):
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.metrics import roc_curve, auc
    if args.use_docking_score:
        X = features
    else:
        X = features[:, :-1]

    y = np.array([1 if x > 0 else 0 for x in label_data])
    n_actives = np.sum(y==1)
    if 0 < args.train_size < 1:
        train_size = args.train_size
    elif args.train_size >= 1:
        train_size = float(args.train_size) / n_actives

    else:
        raise ValueError("train_size must be a positive fraction or active count.")
    if not 0 < train_size < 1:
        raise ValueError("train_size must leave compounds for testing.")
    test_size = 1-train_size

    if args.model == "RF":
        classifier = random_forest(args.random_state, args.nprocs)
    elif args.model == "SVM":
        classifier = support_vector(args.random_state)

    aucs, ef10s, ef1s = [], [], []

    # train-test split, stratified, default:
    # n_splits=5, random_state=None
    splitter = StratifiedShuffleSplit(n_splits=args.n_splits,
                                      random_state=args.random_state,
                                      train_size=train_size,
                                      test_size=test_size)
    for train, test in splitter.split(X, y):
        classifier = classifier.fit(X[train], y[train])
        probas_ = classifier.predict_proba(X[test])
        if args.reverse == True:
            probability = probas_[:, 0]
        else:
            probability = probas_[:, 1]

        fpr, tpr, thresholds = roc_curve(y[test], probability)
        aucs.append(auc(fpr, tpr))

        from calc_ef import calc_ef
        sorted_y = y[test][np.argsort(-probability, kind="stable")]

        ef10 = calc_ef(sorted_y, None if args.active is None else args.active*test_size,
                        None if args.decoy is None else args.decoy*test_size, threshold=0.1)
        ef10s.append(ef10)
        ef1 = calc_ef(sorted_y, None if args.active is None else args.active*test_size,
                        None if args.decoy is None else args.decoy*test_size, threshold=0.01)
        ef1s.append(ef1)

    result = [[np.mean(aucs), np.std(aucs)],
              [np.mean(ef10s), np.std(ef10s)],
              [np.mean(ef1s), np.std(ef1s)]]
    result = pd.DataFrame(result, index=["AUC", "EF10%", "EF1%"],
                          columns=["mean", "std"])
    result.to_csv(args.output, sep=",")

    logger.info("SIEVE: datasize-analysis is finished.")


def cv_accuracy_plot(clf, features, labels, cpd_names, model_name, args):
    from numpy import interp
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_curve, auc

    docking_score = features[:, -1]

    if args.use_docking_score:
        X = features
    else:
        X = features[:, :-1]
    y = labels

    """ k-fold cv, make ROC for each classifier """
    cvs = StratifiedKFold(n_splits=cv, shuffle=True, random_state=args.random_state)

    scores = []
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    aucs = []
    efs = []
    mean_importances = np.array([0.0 for _ in range(X.shape[1])])
    i = 0
    for train, test in cvs.split(X, y):

        clf = clf.fit(X[train], y[train])
        probas_ = clf.predict_proba(X[test])
        cpd_names_ = cpd_names[test]


        # Record the requested class probability.
        if args.reverse == True:
            probas = probas_[:, 0]
        else:
            probas = probas_[:, 1]

        scores.append(pd.DataFrame({"name": cpd_names_, "score": probas}))

        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test], probas)
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, label='%s fold %d (AUC = %0.3f)'
                 % (model_name, i, roc_auc))

        sorted_y = y[test][np.argsort(-probas, kind="stable")]

        from calc_ef import calc_ef
        ef10 = calc_ef(sorted_y, args.active, args.decoy, threshold=0.1)
        ef1 = calc_ef(sorted_y, args.active, args.decoy, threshold=0.01)
        efs.append([ef10, ef1])
        i += 1
        # get feature importance
        if "RF" in model_name:
            try:
                mean_importances += clf.feature_importances_
            except AttributeError:
                import traceback
                traceback.print_exc()
                logger.debug("feature importance is not available. Not Forests?")

    mean_tpr /= cv
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot(mean_fpr, mean_tpr, 'k--',
             label='Mean (AUC = %0.2f)' % mean_auc, lw=1)

    # save
    sorted_scores = pd.concat(scores, ignore_index=True).sort_values(
        "score", ascending=False, kind="stable")
    sorted_scores.to_csv("result_" + model_name + "_scores.csv", header=False, index=False)

    # feature importance for RF
    if "RF" in model_name:
        mean_importances /= float(cv)
        outfile = splitext(args.output)[0] + ".importance"
        np.savetxt(outfile, mean_importances, delimiter=",")

    # glide score, reverse order
    fpr, tpr, thresholds = roc_curve(y, docking_score * (-1))
    docking_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, 'r--', lw=1, label='Glide SP (AUC = %0.3f)' % docking_auc)

    plt.plot([0, 1], [0, 1], '--', lw=1, color=(0.6, 0.6, 0.6), label='Random')

    # formatting
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC: ' + args.title)
    plt.legend(loc="lower right", fontsize=9)
    outfile = splitext(args.output)[0] + "_auc.png"
    plt.savefig(outfile)
    plt.clf()

    aucs.append(mean_auc)
    aucs.append(docking_auc)
    out_auc = np.array(aucs).T
    outfile = splitext(args.output)[0] + "_auc.csv"
    np.savetxt(outfile, out_auc, delimiter=",", fmt="%.3f")

    np_efs = np.array(efs)
    mean_ef10 = np.mean(np_efs[:, 0])
    mean_ef1 = np.mean(np_efs[:, 1])
    efs.append([mean_ef10, mean_ef1])

    docking_sorted_y = y[np.argsort(docking_score, kind="stable")]
    docking_ef10 = calc_ef(docking_sorted_y, args.active, args.decoy, threshold=0.1)
    docking_ef1 = calc_ef(docking_sorted_y, args.active, args.decoy, threshold=0.01)
    efs.append([docking_ef10, docking_ef1])
    out_ef = np.array(efs)
    outfile = splitext(args.output)[0] + "_ef.csv"
    np.savetxt(outfile, out_ef, delimiter=",", fmt="%.3f")

    return aucs, out_ef


def scoring_param_search(title, label_data, interaction_name, features, args):
    if args.zeroneg:
        y = np.array([1 if x > 0 else 0 for x in label_data])
    else:
        y = label_data

    if args.use_docking_score:
        X = features
    else:
        X = features[:, :-1]


    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=cv)

    if args.model == "RF":
        # Random Forest, Grid search by CV
        model = random_forest(args.random_state)
        max_f = min(31, len(X[0, :]))
        param_grid = [{'max_features': list(range(2, max_f)) + ['sqrt']}]

    elif args.model == "SVM":
        model = support_vector(args.random_state)
        param_grid = [{"C": [1, 10, 100], "gamma": [0.01, 0.1, 1]}]

    from sklearn.model_selection import GridSearchCV
    clf = GridSearchCV(model, param_grid, cv=skf,
                                   scoring='roc_auc', n_jobs=args.nprocs)
    clf.fit(X, y)

    # evaluate scores
    results = clf.cv_results_
    grid_scores_df = pd.DataFrame(results["params"])
    grid_scores_df["mean_score_CV"] = results["mean_test_score"]
    grid_scores_df["std_CV"] = results["std_test_score"]

    # output
    grid_scores_df.to_csv(args.output, sep=",", index=False)

    print(clf.best_params_, clf.best_score_)
    score = clf.best_estimator_.predict_proba(X)[:, 1]

    rank = np.argsort(-score, kind="stable")[:args.propose]
    cpd_name = title[rank]
    score = score[rank]
    label = y[rank]


    logger.info('Saved SIEVE-Score.')

    return cpd_name, score, label


def scoring_eval(cpd_names, label_data, interaction_name, features, args):

    if args.zeroneg:
        labels = np.array([1 if x > 0 else 0 for x in label_data])
    else:
        # TODO
        logger.info("not zeroneg is not implemented here.")
        raise NotImplementedError()

    if args.model == "RF":
        classifier = random_forest(args.random_state, args.nprocs)
        model_name = 'SIEVE-Score_RF'

    elif args.model == "SVM":
        classifier = support_vector(args.random_state)
        model_name = 'SIEVE-SVM'

    mean_auc, efs = cv_accuracy_plot(classifier, features, labels, cpd_names, model_name, args)

    if args.model == "RF":
        importance_file = splitext(args.output)[0] + ".importance"
        importance = pd.read_csv(importance_file, header=None, dtype="float64")
        if not args.use_docking_score:
            interaction_name = interaction_name[:-1]
        interaction_name = pd.DataFrame(interaction_name)
        importance = pd.concat([interaction_name, importance], axis=1)
        importance.to_csv(importance_file, sep=",", header=False, index=False)

    return


def scoring_compareSVMRF(title, label_data, interaction_name, features, args):
    if args.zeroneg:
        labels = np.array([1 if x > 0 else 0 for x in label_data])
    else:
        # TODO
        logger.info("not zeroneg is not inplemented here.")
        raise NotImplementedError()

    classifiers = [support_vector(args.random_state),
                   random_forest(args.random_state, args.nprocs)]
    names = ['SIEVE-SVM', 'SIEVE-Score_RF']
    for clf, name in zip(classifiers, names):
        mean_auc, efs = cv_accuracy_plot(clf, features, labels, title, name, args)
        print(mean_auc, efs)
    logger.info('Finished to compare')
    return
