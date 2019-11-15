import numpy as np
import pandas as pd
import logging
from os.path import splitext

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)
cv = 5  # define k of k-fold CV


def lessdata_auc_ef(cpd_names, label_data, interaction_name, features, args):
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.metrics import roc_curve, auc
    from scipy import interp
    
    if args.use_docking_score:
        X = features
    else:
        X = features[:, :-1]

    y = np.array([1 if x > 0 else 0 for x in label_data])
    n_actives = np.sum(y==1)
    if 0 < args.train_size < 1:
        train_size = args.train_size
    elif args.train_size > 1:
        train_size = float(args.train_size) / n_actives

    print(args.train_size, n_actives, train_size)
    if train_size >= 1: 
        quit()
    test_size = 1-train_size

    if args.model == "RF":
        from sklearn.ensemble import RandomForestClassifier as RFC
        classifier = RFC(n_estimators=1000, criterion='gini', max_features=6,
                         random_state=args.random_state)
    elif args.model == "SVM":
        from sklearn.svm import SVC
        classifier = SVC(C=10, kernel="rbf", degree=3, gamma=0.1,
                          cache_size=1000, probability=True)

    model_names = 'SIEVE-Score ' + args.model
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
        cpd_names_ = cpd_names[test]
        if args.reverse == True:
            probability = probas_[:, 0]
        else:
            probability = probas_[:, 1]
        
        score = np.vstack((cpd_names_, probability)).T
        fpr, tpr, thresholds = roc_curve(y[test], probability)
        fpr = [0.0] + fpr
        tpr = [0.0] + tpr
        aucs.append(auc(fpr, tpr))

        from calc_ef import calc_ef
        sorted_y = y[test][probability.argsort()][::-1]
        sorted_cpdnames_ = cpd_names_[probability.argsort()][::-1]
        sorted_probas_ = probas_[probability.argsort()][::-1]
        
        n_active = np.sum(y[test] == 1)
        n_decoy = np.sum(y[test] == 0)
        ef10 = calc_ef(sorted_y, args.active*test_size, args.decoy*test_size, threshold=0.1)
        ef10s.append(ef10)
        ef1 = calc_ef(sorted_y, args.active*test_size, args.decoy*test_size, threshold=0.01)
        ef1s.append(ef1)

    result = [[np.mean(aucs), np.std(aucs)],
              [np.mean(ef10s), np.std(ef10s)],
              [np.mean(ef1s), np.std(ef1s)]]
    result = pd.DataFrame(result, index=["AUC", "EF10%", "EF1%"],
                          columns=["mean", "std"])
    result.to_csv(args.output, sep=",")

    logger.info("SIEVE: datasize-analysis is finished.")


def cv_accuracy_plot(clf, features, labels, cpd_names, model_name, args):
    from scipy import interp
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_curve, auc

    docking_score = features[:, -1]
    print(docking_score)
    if args.use_docking_score:
        X = features
    else:
        X = features[:, :-1]
    y = labels

    """ k-fold cv, make ROC for each classifier """
    cvs = StratifiedKFold(n_splits=cv, shuffle=True)

    scores = np.array([])
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    aucs = []
    efs = []
    mean_importances = np.array([0.0 for _ in range(X.shape[1])])
    i = 0
    for train, test in cvs.split(X, y):
        print(y[train], y[test])
        clf = clf.fit(X[train], y[train])
        probas_ = clf.predict_proba(X[test])
        cpd_names_ = cpd_names[test]
        print(probas_)

        # record scores, having bug for some target?
        if args.reverse == True:
            probas = probas_[:, 0]
        else:
            probas = probas_[:, 1]
            
        score = np.vstack((cpd_names_, probas)).T
        if np.size(scores) > 0:
            scores = np.concatenate([scores, score], axis=0)
        else:
            scores = score

        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test], probas)
        tpr = [0.0] + tpr
        fpr = [0.0] + fpr
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, label='%s fold %d (AUC = %0.3f)'
                 % (model_name, i, roc_auc))

        sorted_score = np.sort(probas)[::-1]
        sorted_y = y[test][probas.argsort()][::-1]

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
    sorted_scores = scores[np.argsort(scores[:, 1])][::-1]
    np.savetxt("result_" + model_name + "_scores.csv", sorted_scores, fmt="%s", delimiter=",")
    
    # feature importance for RF
    if "RF" in model_name:
        mean_importances /= float(cv)
        outfile = splitext(args.output)[0] + ".importance"
        np.savetxt(outfile, mean_importances, delimiter=",")

    # glide score, reverse order
    fpr, tpr, thresholds = roc_curve(y, docking_score * (-1))
    tpr = [0.0] + tpr
    fpr = [0.0] + fpr
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
    
    docking_sorted_y = y[docking_score.argsort()]
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

    docking_score = features[:, -1]
    if args.use_docking_score:
        X = features
    else:
        X = features[:, :-1]


    from sklearn.cross_validation import StratifiedKFold
    skf = StratifiedKFold(y, cv)

    if args.model == "RF":
        # Random Forest, Grid search by CV
        from sklearn.ensemble import RandomForestClassifier
        model = RandomForestClassifier(n_estimators=1000, criterion='gini')
        max_f = min(31, len(X[0, :]))
        param_grid = [{'max_features': list(range(2, max_f)) + ['auto']}]

    elif args.model == "SVM":
        from sklearn.svm import SVC
        model = SVC(C=10, kernel="rbf", degree=3, gamma=0.1, cache_size=1000)
        param_grid = [{""}]

    from sklearn import grid_search
    clf = grid_search.GridSearchCV(model, param_grid, cv=skf,
                                   scoring='roc_auc', n_jobs=args.nprocs)
    clf.fit(X, y)

    # evaluate scores
    grid_scores = [["max_features", "mean_score_CV", "std_CV"]]
    for params, mean_score, all_scores in clf.grid_scores_:
        grid_scores.append([params['max_features'],
                            mean_score, all_scores.std()])
    grid_scores_df = pd.DataFrame(grid_scores[1:],
                                  columns=grid_scores[0])
    print(grid_scores_df)

    # output
    grid_scores_df.to_csv(args.output, delimiter=",")

    print(clf.best_params_, clf.best_score_)
    score = clf.best_estimator_.predict_proba(X)[:, 1]

    rank = np.argsort(score)[::-1][:args.propose]
    cpd_name = title[rank]
    score = score[rank]
    label = y[rank]

    result = np.array(zip(cpd_name, score, label))
    # print(result)
    result = np.dstack((cpd_name, score, label))

    # test
    # np.savetxt(outfile, result, fmt="%s", delimiter=",")
    logger.info('Saved SIEVE-Score.')

    return cpd_name, score, label


def scoring_eval(cpd_names, label_data, interaction_name, features, args):

    if args.zeroneg:
         labels = np.array([1 if x > 0 else 0 for x in label_data])
    else:
        # TODO
        logger.info("not zeroneg is not implemented here.")
        raise NotImplementedError()

    random_state = 0
    if args.model == "RF":
        from sklearn.ensemble import RandomForestClassifier as RFC
        classifier = RFC(n_estimators=1000, criterion='gini', max_features=6, 
                         n_jobs=args.nprocs)
        model_name = 'SIEVE-Score_RF'

    elif args.model == "SVM":
        from sklearn.svm import SVC
        classifier = SVC(C=10, kernel="rbf", degree=3, gamma=0.1,
                         cache_size=1000, probability=True)
        model_name = 'SIEVE-SVM'
    
    mean_auc, efs = cv_accuracy_plot(classifier, features, labels, cpd_names, model_name, args)

    if args.model == "RF":
        importance_file = splitext(args.output)[0] + ".importance"
        importance = pd.read_csv(importance_file, header=None, dtype="float64")
        if not args.use_docking_score:
            interaction_name = interaction_name[:-1]
        interaction_name = pd.DataFrame(interaction_name)
        importance = pd.concat([interaction_name, importance],axis=1, join_axes=[importance.index])
        importance.to_csv(importance_file, sep=",", header=False, index=False)

    return


def scoring_compareSVMRF(title, label_data, interaction_name, features, args):
    outfile = args.output

    if args.zeroneg:
        labels = np.array([1 if x > 0 else 0 for x in label_data])
    else:
        # TODO
        logger.info("not zeroneg is not inplemented here.")
        raise NotImplementedError()

    from sklearn.ensemble import RandomForestClassifier as RFC
    from sklearn.svm import SVC

    classifiers = [SVC(C=10, kernel="rbf", degree=3, gamma=0.1, cache_size=1000, probability=True),
                   RFC(n_estimators=1000, criterion='gini', max_features=6)]
    names = ['SVM', 'RF']
    for clf, name in zip(classifiers, names):
        mean_auc, efs = cv_accuracy_plot(clf, features, labels, name, args)
        print(mean_auc, efs)
    logger.info('Finished to compare')
    quit()
